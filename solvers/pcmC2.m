function [betaPCMMat, sigmaPCMMat,tauPCMMat,out] = pcmC2(Xmat, Yvec, inopts)
% Perspective M-estimators (Combettes and Mueller, 2018,2019)
%
% Generic proximal solver for penalized perspective M-estimators
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: PCM solutions: betaPCMMat, sigmaPCMMat,tauPCMMAT,out

[n, p] = size(Xmat);

% Options for concomitant data fitting term
defopts.objFun = 'Huber';           % Lq, Huber, BerHu, Vapnik
defopts.rho1 = 1.345;               % Scale parameter for Huber/Vapnik/Scaled function
defopts.qPower1 = 2;                % Exponent of the perspective function
defopts.fitLin = 1/2*ones(n,1);     % Factor in front of linear term for concomitant estimation
defopts.fitLB = 0;                  % Lower bound on concomitant variables
defopts.nGroupSiz = n;              % Homoscedastic model has one group of size n

% Options for concomitant penalty fitting term
defopts.penFun = '';                % L2, Lp, Huber/Lp, Vapnik, BerHu
defopts.rho2 = 1.345;               % Scale parameter for Huber/Vapnik/Scaled function
defopts.qPower2 = 2;                % Exponent of the perspective function
defopts.penLin = 1/2*ones(p,1);     % Factor in front of linear term for concomitant estimation
defopts.penLB = 0;                  % Lower bound on concomitant variables
defopts.Lmat = eye(p,p);            % Linear operator on beta vector
defopts.rVec = zeros(p,1);          % Shift in penalization argument
defopts.pGroupSiz = p;              % Homoscedastic penalty has one group of size p

% Options for simple regularization term on beta vector
defopts.regFun = '';%'L1';          % '', L1, L2, L1s (restricted to sum(b)=0 for log-contrast model)
defopts.regWeights = ones(p,1);     % element-wise  weights for the norms (e.g., for the adaptive Lasso)
defopts.Ceq = ones(1,p);            % Matrix needed for log-contrast model constraints
defopts.rhsvec = 0;                 % RHS of the linear equality 

% Options for regularization path
defopts.lenPath = 200;              % Length of lambda path
defopts.delta = 2;                  % Constant for geometric sequence spacing of the regularization path
defopts.activeSet = 1:p;            % Active set of variables
defopts.lamPath = [];               % User input lambda values
defopts.warmstart = 1;              % Option for warmstart over the lambda path

% Options for Douglas Rachford scheme
defopts.abstol = 1e-4;              % DR tolerance (in location), 10*abstol is used for the scale
defopts.maxit = p*1e3;              % p*1e2 number of iteration in DR
defopts.minit = 50;                 % minimum number of iteration in DR (used to prevent premature convergence)
defopts.dr_mu = 1.9;                % mu in DR
defopts.gamma = 1;                  % gamma in DR

% Options for algorithm output
defopts.plotting = 0;               % Plotting of optimization trace
defopts.verbose = 0;

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

if ~isempty(opts.lamPath)
    opts.lenPath = length(opts.lamPath);
    lamVec = opts.lamPath; % Set from user input
    lamVec = sort(lamVec,'descend'); % Sort from largest to smallest IMPORTANT!!!!
else
    lamVec = [];
end

% Length of the regularization path after initialization
nLams = opts.lenPath;

% Geometric sequence for lambda path if not set externally
% Scaling only works for Lasso so far...
if nLams>1 && isempty(lamVec)
    
    % Data-dependent lambda upper bound
    lam_max = n*max(abs(Xmat'*Yvec))./(sqrt(sum(Yvec.^2))*sqrt(n-1));
    
    % Geometric sequence for lambda path
    delta = opts.delta;
    lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1));
    
    if strcmp(opts.penFun,'BerHu')
        lamVec = 0.5*lamVec;
    end
    
end

% Add path back to opts structure (used in DR for warmstart checks)
opts.lamPath = lamVec;


% Check grouping of data and predictors
if sum(opts.nGroupSiz)~=n
    error('Grouping information in data fitting misspecified')
end

if sum(opts.pGroupSiz)~=p
    error('Grouping information in penalty misspecified')
end

% Data groups (for heteroscedastic estimation)
nGroups = length(opts.nGroupSiz);

nGroupVec = zeros(n,2);
sInds = 1;
eInds = 0;

% Groups are in linear order
for i=1:nGroups
    eInds = eInds + opts.nGroupSiz(i);
    nGroupVec(sInds:eInds,2) = i;
    sInds = eInds+1;
end

% Weight on linear part
nGroupVec(:,1) = opts.fitLin;

% Regression groups (for structured (sparse) estimation)
pGroups = length(opts.pGroupSiz);

pGroupVec = zeros(p,2);
sInds = 1;
eInds = 0;

for i=1:pGroups
    eInds = eInds + opts.pGroupSiz(i);
    pGroupVec(sInds:eInds,2) = i;
    sInds = eInds+1;
end

% Weight on linear part
pGroupVec(:,1) = opts.penLin;

%% Set up proximity operators

% Set proximity and projection operators from input options
objFun = opts.objFun;
penFun = opts.penFun;
regFun = opts.regFun;

% Scale, exponent, and linear offset for the concomitant data fitting variable
rho1 = opts.rho1;
qPower1 = opts.qPower1;
fitLin = opts.fitLin;
fitLB = opts.fitLB; % lower bound for "smooth" objective

%  Scale, exponent, and linear offset for the predictor variable
rho2 = opts.rho2;
qPower2 = opts.qPower2;
penLin = opts.penLin;
penLB = opts.penLB;

% Additional model for structured sparsity
Lmat = opts.Lmat;
rVec = opts.rVec;

% Element-wise weights for simple penalty (for adaptive L1, ordered L1,...)
regWeights = opts.regWeights;

% Log-contrast model parameters
Ceq = opts.Ceq;           
rhsvec = opts.rhsvec;  

% Scaling parameter for Douglas Rachford
gamma = opts.gamma;

%% Proximity operators for perspective functions

% Data fitting term
if (strcmp(objFun,'Huber') && qPower1==2)
    prox_phi = @(sig_i,b_i,gamma) proxg2Huber(sig_i,b_i,rho1,fitLin(1),gamma);
elseif (strcmp(objFun,'Huber') && qPower1~=2)
    prox_phi = @(sig_i,b_i,gamma) proxgqHuber(sig_i,b_i,rho1,fitLin(1),qPower1,gamma);
elseif (strcmp(objFun,'Lq') && qPower1==2) % Special case scaledLasso
    kappa = 2; % Factor in denominator;
    prox_phi = @(sig_i,b_i,gamma) proxg2Lq(sig_i,b_i,kappa,fitLin(1),gamma);
elseif (strcmp(objFun,'Lq') && qPower1~=2)
    kappa = 2; % Factor in denominator of qth-power
    prox_phi = @(sig_i,b_i,gamma) proxgqLq(sig_i,b_i,kappa,fitLin(1),qPower1,gamma);
elseif (strcmp(objFun,'BerHu') && qPower1==2)
    kappa = 1; % Slope of the L1 part of BerHu; special case standard BerHu
    prox_phi = @(sig_i,b_i,gamma) proxg2BerHu(sig_i,b_i,rho2,fitLin(1),kappa,gamma);
elseif (strcmp(objFun,'BerHu') && qPower1~=2)
    kappa = 1; % Slope of the L1 part of BerHu
    prox_phi = @(sig_i,b_i,gamma) proxgqBerHu(sig_i,b_i,rho2,fitLin(1),kappa,qPower1,gamma);
else
    error([objFun, ' not yet available'])
end

% Penalty term
if (strcmp(penFun,'BerHu') && qPower2==2)
    % Slope of the L1 part
    kappa = 1;
    prox_psi = @(tau_i,b_i,gam_lam) proxg2BerHu(tau_i,b_i,rho2,penLin(1),kappa,gam_lam);
elseif (strcmp(penFun,'BerHu') && qPower2~=2)
    % Slope of the L1 part
    kappa = 1;
    prox_psi = @(tau_i,b_i,gam_lam) proxgqBerHu(tau_i,b_i,rho2,penLin(1),kappa,qPower2,gam_lam);
elseif ((strcmp(penFun,'Lq') || strcmp(penFun,'Organic')) && qPower2==2) % Special case scaledLasso
    kappa = 2; % Factor in denominator;
    prox_psi = @(tau_i,b_i,gam_lam) proxg2Lq(tau_i,b_i,kappa,penLin(1),gam_lam);
elseif ((strcmp(penFun,'Lq') || strcmp(penFun,'Organic')) && qPower2~=2)
    kappa = 2; % Factor in denominator of qth-power
    prox_psi = @(tau_i,b_i,gam_lam) proxgqLq(tau_i,b_i,kappa,penLin(1),qPower2,gam_lam);
elseif (strcmp(penFun,'Lsum') && qPower2==2)
    kappa = 2; % Factor in denominator of qth-power
    prox_psi = @(tau_i,b_i,gam_lam) proxgqLsum(tau_i,b_i,kappa,penLin,qPower2,gam_lam);
elseif (strcmp(penFun,'L1')) 
    prox_psi = @(tau_i,b_i,gam_lam) proxgL1(tau_i,b_i,gam_lam);
elseif (strcmp(penFun,''))
    % Identity mapping
    prox_psi = @(tau_i,b_i,gam_lam) proxgqId(tau_i,b_i,gam_lam);
else
    error([penFun, ' not yet available'])
end

%% Simple proximity operator for Theta function (e.g., L1,L2,adaptive L1, ordered L1, ...)
if (strcmp(regFun,'L1'))
    prox_theta = @(b,lam,gamma) proxL1(b,lam,gamma,regWeights); %regWeights for adaptive Lasso
elseif (strcmp(regFun,'L2'))
    prox_theta = @(b,lam,gamma) proxL2(b,lam,gamma);
elseif (strcmp(regFun,'L1s'))
    prox_theta = @(b,lam,gamma) proxL1s(b,lam,gamma,regWeights);
elseif (strcmp(regFun,'L1o'))
    prox_theta = @(b,lam,gamma) proxL1o(b,lam,gamma,regWeights);
elseif (strcmp(regFun,'Ceq'))
    PCeq = Ceq'*pinv(Ceq*Ceq');
    prox_theta = @(b,lam,gamma) projCeq(b,PCeq,Ceq,rhsvec);
elseif (strcmp(regFun,''))
    prox_theta = @(b,lam,gamma) b;
end

%% Projection operators for concomitant variables
if (strcmp(objFun,'Huber') || strcmp(objFun,'Lq') || strcmp(objFun,'BerHu'))
    projD = @(sVec) projC(sVec,nGroupVec(:,2),fitLB);
else
    projD = @(sVec) sVec;
end

if (strcmp(penFun,'BerHu')) || (strcmp(penFun,'Lq'))
    projE = @(tauVec) projA(tauVec,pGroupVec(:,2),penLB);
elseif (strcmp(penFun,'Organic'))
    projE = @(tauVec) projO(tauVec,pGroupVec(:,2),penLB);
else
    projE = @(tauVec) tauVec;
end


% Create storage for solution vectors
%betaPCMMat = zeros(p,nLams);
%sigmaPCMMat = zeros(nGroups,nLams);
%tauPCMMat = zeros(pGroups,nLams);

% Temporary use of full concomitant vector independent of group identity
betaPCMMat = zeros(p,nLams);
sigmaPCMMat = zeros(n,nLams);
tauPCMMat = zeros(p,nLams);

% Create storage for function values and timings
runTimes = zeros(1, nLams);

% Loop over regularization path
for i = 1:nLams
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i);
    if opts.verbose
        disp(['lambda: ',num2str(lam_i)])
    end
    
    % Adaptive mu schedule for path problems
    % opts.dr_mu = max(0.1,opts.dr_mu*0.95);
    
    % Generalized Douglas-Rachford
    [bVec,sVec,tVec] = dougRach(Yvec,Xmat,rVec,Lmat,projD,projE,prox_theta,prox_phi,prox_psi,lam_i,opts);
    
    betaPCMMat(:,i) = bVec;
    sigmaPCMMat(:,i) = sVec;
    tauPCMMat(:,i) = tVec;
    
    runTimes(1,i) = toc;
    
    % Break if lower bound on lambda is reached
    if sum(abs(bVec)>1e-2)>=n
        remLen = (nLams-i);
        betaPCMMat(:,i+1:end)=repmat(bVec,1,remLen);
        sigmaPCMMat(:,i+1:end)=repmat(sVec,1,remLen);
        tauPCMMat(:,i+1:end)=repmat(tVec,1,remLen);
        runTimes(1,i+1:end) = 0;
        disp('Lower bound on lambda reached')
        break;
    end
    
    %figure(234)
    %plot(lamVec(1:i),betaPCMMat(:,1:i)')
    %grid on
    %drawnow
    
end

out.opts = opts;
out.runTimes = runTimes;
out.X = Xmat;
out.Y = Yvec;
out.lamPath = lamVec;
end

% Douglas-Rachford algorithm
function [bVec,sVec,tVec] = dougRach(Yvec,Xmat,rVec,Lmat,projD,projE,prox_theta,prox_phi,prox_psi,lam,dropts)

% Persistent variables for warm starts
persistent A Q c_bk x_sk x_tk x_bk h_sk h_tk h_bk d_sk d_tk d_bk;

% Dimensionality of the input
[n,p] = size(Xmat);
[n1,p1] = size(Lmat); % Not yet used

% DR parameters
maxit = dropts.maxit;
minit = dropts.minit;
abstol = dropts.abstol;
mu_k = dropts.dr_mu;
warmstart = dropts.warmstart;
gamma = dropts.gamma;

% Upper bound on lambda
lamMax = max(dropts.lamPath);

% Scaling for perspective function of the penalty
%gam_lamMax = gamma*lamMax;

% Rescaling of gamma such that ratio is fixed
%gamma = gam_lamMax/(2*lam);

gam_lam = gamma*lam;


if lam==lamMax || ~warmstart
    % This is the first time we call dougRach...
    
    % Lmat could be used for fusion penalties etc.
    
    % Concatenate both linear transforms (X and L)
    A = [Xmat;Lmat];
    
    % Construct Q matrix
    Q = A'/(eye(n+p,n+p)+A*A');
    
    % All splitting variables
    
    % Decompose c_k
    s_k = ones(n,1);
    t_k = ones(p,1);
    c_bk = zeros(n+p,1);
    
    % Decompose x_k
    x_sk = zeros(n,1);
    x_tk = zeros(p,1);
    x_bk = zeros(p,1);
    
    % Decompose h_k
    h_sk = ones(n,1);
    h_tk = ones(p,1);
    h_bk = zeros(n+p,1);
    
    % Decompose d_k
    d_sk = ones(n,1);
    d_tk = ones(p,1);
    d_bk = zeros(n+p,1);
    
end

b_k_1 = randn(p,1);
s_k_1 = randn(n,1);

% Main Douglas-Rachford loop
for kk=1:maxit
    
    % Additional splits
    q_sk = x_sk - h_sk;
    q_tk = x_tk - h_tk;
    q_bk = A*x_bk - h_bk;
    
    % Convergent sequences
    s_k = x_sk - q_sk/2;
    t_k = x_tk - q_tk/2;
    b_k = x_bk - Q*q_bk;
    
    % Projection operators
    z_sk = projD(2*s_k-x_sk);
    z_tk = projE(2*t_k-x_tk);
    
    % Simple proximity operator on b (e.g., L1) or projection operator in
    % case of convex constraints on b
    z_bk = prox_theta(2*b_k-x_bk,lam,gamma);
    
    % Advance x variables
    x_sk = x_sk + mu_k*(z_sk-s_k);
    x_tk = x_tk + mu_k*(z_tk-t_k);
    x_bk = x_bk + mu_k*(z_bk-b_k);
    
    % Proximity operations for the data fitting term
    for i=1:n
        
        % Update c_b vector
        c_ik = Xmat(i,:)*b_k;
        c_bk(i) = c_ik;
        
        eta_ik = h_sk(i);
        h_ik = h_bk(i);
        
        [prox1,prox2] = prox_phi(2*s_k(i)-eta_ik,(2*c_ik-h_ik)-Yvec(i),gamma);
        delt_d_ik = [0,Yvec(i)] + [prox1,prox2]; %
        d_sk(i) = delt_d_ik(1); % delt_ik in pseudo-code
        d_bk(i) = delt_d_ik(2); % d_ik;
        
    end
    
    % Proximity operations for the penalty term
    for i=1:p
        
        % only a temp variable
        c_Nik = Lmat(i,:)*b_k;
        c_bk(n+i) = c_Nik;
        
        eta_Nik = h_tk(i);
        h_Nik = h_bk(n+i);
        
        [prox1,prox2] = prox_psi(2*t_k(i)-eta_Nik,2*c_Nik-h_Nik-rVec(i),gam_lam);
        delt_d_Nik = [0,rVec(i)] + [prox1,prox2];
        d_tk(i) = delt_d_Nik(1);   % delt_ik in pseudo-code
        d_bk(n+i) = delt_d_Nik(2); % d_ik;
        
    end
    
    h_sk = h_sk + mu_k*(d_sk-s_k);
    h_tk = h_tk + mu_k*(d_tk-t_k);
    h_bk = h_bk + mu_k*(d_bk-c_bk);
    
    % Check convergence of the location iterates
    norm_bk = sqrt(sum((abs(b_k_1)-abs(b_k)).^2));
    
    % Check convergence of the scale iterates
    norm_sk = sqrt(sum((abs(s_k_1)-abs(s_k)).^2));

    if ~mod(kk,5e3) && dropts.verbose
        disp(['Error at iteration: ',num2str(kk),': ',num2str(norm_bk)]);
    end
    
    if (norm_bk<abstol) && (norm_sk<10*abstol) && (kk>minit)
        
        disp(['Convergence to l-tol = ', num2str(norm_bk), ' at iteration: ',num2str(kk)])
        disp(['Convergence to s-tol = ', num2str(norm_sk), ' at iteration: ',num2str(kk)])
        break
    end
    
    deltakk = 1e2;
    
    % Plotting for code checking
    if ~mod(kk,deltakk)
        
%         figure(23);
%         plot(kk,log10(abs(c_bk)),'k.','MarkerSize',20)
%         hold on
%         plot(kk,log10(abs(b_k)),'r.','MarkerSize',20)
%         ylabel('Regression vector')
%         drawnow
        
%         figure(42);
%         plot(kk,s_k,'.','MarkerSize',20)
%         hold on
%         ylabel('Concomitant t')
%         drawnow
%         
%         figure(423);
%         semilogy(kk,norm_bk,'.','MarkerSize',20)
%         hold on
%         ylabel('Norm of b iterate')
%         drawnow
        
    end
    
    % Store previous current b_k and s_k for convergence
    b_k_1 = b_k;
    s_k_1 = s_k;

    
    % Adapting mu
    %mu_k = (1-1e-3)*mu_k;
    
end

% Sparsify the b_k estimate using the splitted counterpart x_bk
%if ~strcmp(dropts.regFun,'Ceq')
    b_k(z_bk(:)==0)=0;
%end

% Final estimate
sVec = s_k;
tVec = t_k;
bVec = b_k;

end



% %%
% function eta_yProx = proxgpVapnik(eta_y,epsV,gamma,Yvec)
% % Perspective of Vapnik in (n+1) dimensions
%
% eta = eta_y(1);
% y = eta_y(2:end)-Yvec;
%
% % Dimension of the problem
% n = length(y);
%
% % Norm of y
% y_norm2 = sum(y.^2);
% y_norm = sqrt(y_norm2);
%
% % Check five cases for prox calculation
%
% % Case 1
% if (eta+epsV*y_norm<=0) && (y_norm <= gamma)
%     etaProx=0;
%     yProx=zeros(n,1);
%     % Case 2
% elseif (eta<=-gamma*epsV) && (y_norm > gamma)
%
%     etaProx = 0;
%     yProx = y-gamma*y./y_norm;
%     % Case 3
% elseif eta > -gamma*epsV && y_norm > epsV*eta + gamma*(1+epsV^2);
%
%     etaProx = eta + gamma*epsV;
%     yProx = y - gamma*y./y_norm;
%
%     % Case 4
% elseif y_norm>-eta/epsV && epsV*eta<= y_norm  && y_norm <=epsV*eta + gamma*(1+epsV^2)
%
%     etaProx = (eta+epsV*y_norm)/(1+epsV^2);
%     yProx = epsV*(eta+epsV*y_norm)*y/(y_norm*(1+epsV^2));
%
%     % Case 5
% elseif eta >= 0 && y_norm <= epsV*eta;
%
%     etaProx = eta;
%     yProx = y;
%
% else
%
%     warning('Case not covered')
%
%     etaProx = eta;
%     yProx = y;
%
% end
%
% eta_yProx = [0;Yvec] + [etaProx;yProx];
%
% end
%
%
% %%
% function [eta_yProx] = proxgsVapnik(eta_y,epsV,gamma,Yvec)
% % Perspective of Vapnik in 2 dimensions
%
% % Number of samples
% n = length(Yvec);
%
% etaVec = eta_y(1:n);
%
% % Shift by data vector
% y = (eta_y(n+1:end)-Yvec);
%
% % Number of samples
% n = length(y);
%
% etaProx = zeros(n,1);
% yProx = zeros(n,1);
%
% % Iterate over all data points
% for i=1:n
%
%     % ith component
%     y_i = y(i);
%
%     eta_i = etaVec(i);
%
%     % Check five cases for prox calculation
%     yabs = abs(y_i);
%
%     if (eta_i+epsV*yabs<=0) && (yabs <= gamma)
%         etaProx(i)=0;
%         yProx(i) = 0;
%
%     elseif (eta_i<=-gamma*epsV) && (yabs > gamma)
%
%         etaProx(i) = 0;
%         yProx(i) = y_i-gamma*sign(y_i);
%
%     elseif eta_i > -gamma*epsV && yabs > epsV*eta_i + gamma*(1+epsV^2);
%
%         etaProx(i) = eta_i + gamma*epsV;
%         yProx(i) = y_i - gamma*sign(y_i);
%
%     elseif yabs>-eta_i/epsV && epsV*eta_i<= yabs  && yabs <= (epsV*eta_i + gamma*(1+epsV^2))
%
%         etaProx(i) = (eta_i+epsV*yabs)/(1+epsV^2);
%         yProx(i) = epsV*(eta_i+epsV*yabs)*y_i/(yabs*(1+epsV^2));
%
%         % Case 5
%     elseif eta_i >= 0 && yabs <= epsV*eta_i;
%
%         etaProx(i) = eta_i;
%         yProx(i) = y_i;
%
%     else
%
%         warning('Case not covered')
%
%         etaProx(i) = eta_i;
%         yProx(i) = y_i;
%
%     end
% end
%
% % % Shift by data vector and return (n + n) dimensional vector
% eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];
%
% end
%
%

% ---------------------------------------------------------------
% FUNCTIONS BELOW ARE TAKEN FROM Niko Hansen's CMA-ES code
% ---------------------------------------------------------------
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%


if nargin < 2 || isempty(defopts) % no default options available
    opts=inopts;
    return;
elseif isempty(inopts) % empty inopts invoke default options
    opts = defopts;
    return;
elseif ~isstruct(defopts) % handle a single option value
    if isempty(inopts)
        opts = defopts;
    elseif ~isstruct(inopts)
        opts = inopts;
    else
        error('Input options are a struct, while default options are not');
    end
    return;
elseif ~isstruct(inopts) % no valid input options
    error('The options need to be a struct or empty');
end

opts = defopts; % start from defopts
% if necessary overwrite opts fields by inopts values
defnames = fieldnames(defopts);
idxmatched = []; % indices of defopts that already matched
for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
        for i = 1:size(defnames, 1)
            idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
        end
    else
        idx = strncmp(lower(defnames), lower(name), length(name));
    end
    if sum(idx) > 1
        error(['option "' name '" is not an unambigous abbreviation. ' ...
            'Use opts=RMFIELD(opts, ''' name, ...
            ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
        defname  = defnames{find(idx)};
        if ismember(find(idx), idxmatched)
            error(['input options match more than ones with "' ...
                defname '". ' ...
                'Use opts=RMFIELD(opts, ''' name, ...
                ''') to remove the field from the struct.']);
        end
        idxmatched = [idxmatched find(idx)];
        val = getfield(inopts, name);
        % next line can replace previous line from MATLAB version 6.5.0 on and in octave
        % val = inopts.(name);
        if isstruct(val) % valid syntax only from version 6.5.0
            opts = setfield(opts, defname, ...
                getoptions(val, getfield(defopts, defname)));
        elseif isstruct(getfield(defopts, defname))
            % next three lines can replace previous three lines from MATLAB
            % version 6.5.0 on
            %   opts.(defname) = ...
            %      getoptions(val, defopts.(defname));
            % elseif isstruct(defopts.(defname))
            warning(['option "' name '" disregarded (must be struct)']);
        elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
            opts = setfield(opts, defnames{find(idx)}, val);
            % next line can replace previous line from MATLAB version 6.5.0 on
            % opts.(defname) = inopts.(name);
        end
    else
        warning(['option "' name '" disregarded (unknown field name)']);
    end
end
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function res = isoctave
% any hack to find out whether we are running octave
s = version;
res = 0;
if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
end
end
