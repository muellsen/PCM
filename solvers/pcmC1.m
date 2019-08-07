function [betaPCMMat, sigmaPCMMat,funPCMMat,out] = pcmC1(Xmat, Yvec, inopts)
% Penalized Concomitant M-estimators (Combettes and Mueller, 2017)
% with more concomitant variables
% Generic proximal solver for penalized concomitant estimators
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: PCM solutions: betaPCMMat, sigmaPCMMat,funPCMMat,out

[n, p] = size(Xmat);

% Options for concomitant data fitting term
defopts.objFun = 'Huber';           % L2, Lp, Huber/Lp, Vapnik, BerHu
defopts.rho1 = 1.345;               % Scale parameter for Huber/Vapnik/Scaled function
defopts.fitLin = 1/2*ones(n,1);     % factor in front of linear term for concomitant estimation
defopts.nGroupSiz = n;              % Homoscedastic model has one group of size n
defopts.qPower1= 3/2;               % qth power for the Lq data fitting

% Options for concomitant penalty term
defopts.penFun = 'L1';              % L1, L2, Lp, Huber/Lp, BerHu
defopts.rho2 = 1.345;               % Scale parameter for Huber/Vapnik/Scaled function
defopts.penLin = 1/2*ones(p,1);     % factor in front of linear term for concomitant estimation
defopts.LinOperator = eye(p,p);     % Linear operator on beta vector
defopts.pGroupSiz = p;              % Homoscedastic penalty has one group of size p
defopts.qPower2= 3/2;               % qth power for the Lq penalty (aka Bridge)

% Options for regularization path
defopts.lenPath = 200;              % Length of lambda path
defopts.delta = 2;                  % Constant for geometric sequence spacing of the regularization path
defopts.activeSet = 1:p;            % Active set of variables
defopts.lambda = [];                % Single lambda value

% Options for proximal scheme
defopts.abstol = 1e-4;              % DR tolerance
defopts.maxit = p*5e2;              % p*1e3 number of iteration in DR
defopts.dr_mu = 1.9;                % mu in DR
defopts.gamma = n/p;                  % gamma in DR

% Options for algorithm output
defopts.plotting = 0;               % Plotting of trace
defopts.verbose = 0;


% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

if ~isempty(opts.lambda)
    opts.lenPath = length(opts.lambda);
    lamVec = opts.lambda;
else
    lamVec = [];
end


% Length of the regularization path
nLams = opts.lenPath;

% Check grouping of data
if sum(opts.nGroupSiz)~=n
    error('Grouping information misspecified')
end

nGroups = length(opts.nGroupSiz);
pGroups = length(opts.pGroupSiz);

groupVec = zeros(n,2);
sInds = 1;
eInds = 0;
for i=1:nGroups
    eInds = eInds + opts.nGroupSiz(i);
    groupVec(sInds:eInds,2) = i;
    sInds = eInds+1;
end

% Weight on linear part
groupVec(:,1) = opts.fitLin;

% Create storage for 2 nA beta vectors, function values and timings
betaPCMMat = zeros(p,nLams);
sigmaPCMMat = zeros(nGroups,nLams);
tauPCMMat = zeros(pGroups,nLams);
funPCMMat = zeros(1, nLams);
runTimes = zeros(1, nLams);
betaOrbitCell = cell(1,nLams);

% Algorithm stopping criteria
maxit = opts.maxit;
abstol = opts.abstol;

% DR gamma
gamma = opts.gamma;

% Concomitant scales
rho1 = opts.rho1;

rho2 = opts.rho2;

% Linear coefficients + grouping

% Data-dependent lambda upper bound
%lam_max = max(Xmat'*Yvec)./(sqrt(sum(Yvec.^2))*sqrt(n));
%lam_max = max(Xmat'*Yvec)./(sqrt(sum(Yvec.^2)));

% For Huber penalty
rhoInds = (abs(Yvec)>rho1);
nrhoInds = (abs(Yvec)<=rho1);
lam_max = n*max(Xmat'*Yvec)./sum([abs(Yvec(rhoInds));(Yvec(nrhoInds)).^2/2]);

%lam_max = p*lam_max

% Initial beta (solution vector)
beta0 = zeros(p,1);

% Initial lambda
lam_0 = lam_max;

% Geometric sequence for lambda path if not set externally
if nLams>1 && isempty(lamVec)
    delta = opts.delta;
    lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1));
end

% Initial sigma (concomitant variable)
%sigma0 = sqrt(sum(Yvec.^2))/sqrt(n);
sigma0 = 0.5;

% Inital tau
tau0 = 0.5;

% % Precompute properties of the data ONCE
% XY = Xmat'*Yvec;
% XX = Xmat'*Xmat;
% sum_X2 = sum(Xmat.^2);

% Loop over regularization path
for i = 1:nLams
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i);
    
    disp(['lambda: ',num2str(lam_i)])
    
    % Douglas-Rachford
    [beta_K,beta_KDual,betaOrbit,M,R,flag] = dougRach(Yvec,Xmat,beta0,groupVec,sigma0,tau0,rho1,rho2,gamma,lam_i,opts);
    
    % DR returned feasible values
    if flag==0
        
        betaPCMMat(:,i) = beta_K((nGroups+1):end,1);
        sigmaPCMMat(1:nGroups,i) = beta_K(1:nGroups);
        
        if p<100
            betaOrbitCell{1,i} = betaOrbit;
        end
        runTimes(1,i) = toc;
        
        % Exit when Vapnik is violated
        %     if norm(beta_K(2:end,1))<1e-1
        %         disp(['Upper bound reached at $\lambda$: ',num2str(lam_i)])
        %         break
        %     end
        
        %Warm start
        sigma0 = mean(beta_K(1:nGroups));
        beta0 = beta_K((nGroups+1):end,1);
        
        % Something in DR went wrong
    elseif flag==1 || flag==2
        
        for k=i:nLams
            betaPCMMat(:,k) = betaPCMMat(:,i-1);
            sigmaPCMMat(1:nGroups,k) = sigmaPCMMat(1:nGroups,i-1);
        end
        runTimes(1,i:end) = 0;
        break;
    end
    
end

out.opts = opts;
out.runTimes = runTimes;
out.X = Xmat;
out.Y = Yvec;
out.lamPath = lamVec;
out.betaOrbitCell = betaOrbitCell;
out.M = M;
out.R = R;
end

% Douglas-Rachford algorithm
function [betaPCM,betaPCMDual,betaOrbit,M,R,flag] = dougRach(Yvec,Xmat,beta0,groupVec,sigma0,tau0,rho1,rho2,gamma,lam,dropts)

maxit = dropts.maxit;
abstol = dropts.abstol;
mu_k = dropts.dr_mu;

% In case Lq norm is used
qPower1 =dropts.qPower1;
qPower2 =dropts.qPower2;

[n,p] = size(Xmat);

nGroups = max(groupVec(:,2));
groupSInds = zeros(nGroups,1);
for i = 1:nGroups
    groupSInds(i) = find(groupVec(:,2)==i,1);
end

% n+n x n+p matrix
M = zeros(n+n,p+n);
M(1:n,1:n) = eye(n);
M(n+1:n+n,n+1:n+p) = Xmat;

% Matrices
Id_n = eye(n+n);
R = M'/(Id_n+M*M');

% Initialize DR sequences

y0 = Yvec;

sigma0n = sigma0*ones(n,1);

x_k = [sigma0n;beta0]; % n+p vector
y_k = [sigma0n;y0];    % n+n vector

% Help variable
z_k1 = [sigma0n;beta0]; % n+p vector
t_k = zeros(2*n,1);

betaOrbit = zeros(p+n,maxit);

K_old = 1;
w_old = x_k;

fitstr = dropts.objFun;
penstr = dropts.penFun;


% Main DR loop
for kk=1:maxit
    
    q_k = M*x_k - y_k;
    w_k = x_k - R*q_k;
    r_k = M*w_k;
    
    % Projection for groups
    %if strcmp(fitstr,'TREX')
    %
    %    z_k(1:n,1) = projTREX(2*w_k(1:n)-x_k(1:n),2*w_k(n+1:n+p)-x_k(n+1:n+p),groupVec(:,2),groupVec(:,1),gamma,Xmat,Yvec);
    %
    %else
    z_k(1:n,1) = projD(2*w_k(1:n)-x_k(1:n),groupVec(:,2),groupVec(:,1),gamma);
    %end
    
    if strcmp(penstr,'L1')
        
        % Soft thresholding
        z_k(n+1:n+p,1)  = proxST(2*w_k(n+1:n+p)-x_k(n+1:n+p),lam,gamma);
        
    elseif strcmp(penstr,'L1c')
        
        % Soft thresholding
        z_k(n+1:n+p,1)  = proxSTc(2*w_k(n+1:n+p)-x_k(n+1:n+p),lam,gamma);
        
    elseif strcmp(penstr,'L2')
        
        % Ridge regularization
        z_k(n+1:n+p,1) = proxL2(2*w_k(n+1:n+p)-x_k(n+1:n+p),lam,gamma);
        
    end
    
    if strcmp(fitstr,'Huber')
        
        % Data fitting
        t_k = proxgsHuber(2*r_k-y_k,rho1,gamma,Yvec);
        
    elseif strcmp(fitstr,'Vapnik')
        t_k = proxgsVapnik(2*r_k-y_k,rho1,gamma,Yvec);
        
    elseif strcmp(fitstr,'BerHu')
        %t_k = proxgsBerHu(2*r_k-y_k,rho1,gamma,Yvec);
        alpha = 0;
        kappa = 1;
        t_k = proxgqsBerHu(2*r_k-y_k,rho1,alpha,kappa,qPower1,gamma,Yvec);
        %sum(abs((t_k - t_ktemp)))
        %pause
    elseif strcmp(fitstr,'L2') %Scaled Lasso
        
        alpha = 2/n;
        
        temp = 2*r_k-y_k;
        eta_y = temp(n:end,1);
        
        temp2 = proxgpL2(eta_y,alpha,gamma,Yvec);
        
        t_k(n:end,1) = temp2;
        t_k(1:n-1) = t_k(n);
        
    elseif strcmp(fitstr,'Lq')
        
        alpha = 2/n^(qPower1/2);
        
        temp = 2*r_k-y_k;
        
        
        eta_y = temp(n:end,1);
        
        %temp2 = proxgLq(eta_y,alpha,qPower1,gamma,Yvec);
        
        %t_k(n:end,1) = temp2;
        %t_k(1:n-1) = t_k(n);
        
        %temp=(x_k)';
        %temp(1:n)
        %pause
        t_k = proxgsLq(2*r_k-y_k,alpha,qPower1,gamma,Yvec);
        
    else
        error('Fitting unknown')
        
    end
    
    x_k1 = x_k + mu_k*(z_k-w_k);
    y_k1 = y_k + mu_k*(t_k-r_k);
    
    
    xx_k = [w_k;r_k];
    zz_k = [z_k;t_k];
    
    % Check convergence of iterates
    norm_xz = norm(xx_k-zz_k);
    
    if ~mod(kk,5e3)
        disp(['Error at iteration: ',num2str(kk),': ',num2str(norm_xz)]);
    end
    
    if (norm_xz<abstol) && (kk>5e1)
        x_k = x_k1;
        y_k = y_k1;
        
        disp(['Early convergence at iteration: ',num2str(kk)])
        break
    end
    
    z_k1 = z_k;
    
    betaOrbit(:,kk) = z_k;
    x_k = x_k1;
    y_k = y_k1;
    
    % Adapting mu
    %mu_k = (1-1e-3)*mu_k;
    
    if dropts.plotting
        
        if ~mod(kk,2e1)
            semilogy([K_old,kk],abs([w_old(n+1:end),w_k(n+1:end)]),'k-','LineWidth',5)
            grid on
            hold on
            drawnow
            w_old = w_k;
            K_old = kk;
        end
    end
end


if (z_k(1)~=z_k(2))
    z_k(1)
    z_k(2)
    warning('Averaging in concomitant variable did not work')
end

% Final estimate (sigma is first nGroups component, beta follows)

betaPCM = [z_k(groupSInds);z_k(n+1:end)];
betaPCMDual = [y_k(groupSInds);y_k(n+1:end)];

%disp(['First concomitant: ',num2str(z_k(1))])

if z_k(1)<0
    warning('Concomitant variable is negative')
    flag = 1;
elseif abs(z_k(1))<1e-12
    warning('Concomitant variable is zero')
    flag = 2;
else
    flag = 0;
end

end

function sProj = projD(sVec,gV,shV,gam)
% Project and shift for (group) concomitant variables

n = size(gV,1);
nG = max(gV);

sProj = sVec;

% Projection for the n concomitant variables with linear part shV
for i=1:nG
    currInds = (gV==i);
    sProj(currInds) = mean(sVec(currInds))-gam*shV(currInds,1);%ones(n,1);%*ones(n,1)-a;
end
end

function sProj = projTREX(sVec,b,gV,shV,gam,X,Y)
% Project and shift for (group) concomitant variables
% in the TREX model

n = size(gV,1);
nG = max(gV);

sProj = sVec;

tInd = 1;

% Projection for the n concomitant variables with linear part shV
for i=1:nG
    currInds = find(gV==i);
    currS = max(X(currInds,tInd)'*(X(currInds,:)*b-Y(currInds)),0);
    sProj(currInds) = currS*ones(length(currInds),1)-gam*shV(currInds,1);%ones(n,1);%*ones(n,1)-a;
end
end


%%
function betaT = proxST(b,lam,gam)
% Soft thresholding operation on the n+1:n+p components

lam_gam = lam*gam;

% Soft thresholding of the predictors
betaT = sign(b).*max(abs(b) - lam_gam,0);

end

%%
function betaT = proxSTc(b,lam,gam)
% Soft thresholding operation on the n+1:n+p components

lam_gam = lam*gam;

% Soft thresholding of the predictors and zero sum constraint
betaT = prox_l1sum(b,lam_gam,0);

end


%%
function betaT = proxL2(b,lam,gam)
% Proximity operator for L2

lam_gam = lam*gam;

betaT = (1-lam_gam./norm(b))*b;

end

%%
function eta_yProx = proxL2_dep(eta_y,y_hat,gamma)
% Proximity operator for L2
etaProx = eta_y(1);
y = eta_y(2:end);


% First variant
%y = (y-y_hat);
%yProx = y/(1+gamma);

% Second variant
y_hat = y_hat*gamma;
yProx = y/(1+gamma)+y_hat/(1+gamma);

eta_yProx = [etaProx;yProx];

end

%%
function [eta_yProx,caseInds] = proxgpHuber(eta_y,rho,gamma,Yvec)
% Perspective of Huber in (n+1) dimensions

eta = eta_y(1);
y = eta_y(2:end)-Yvec;

% Dimension of the problem
n = length(y);

% Norm of y
y_norm2 = sum(y.^2);
y_norm = sqrt(y_norm2);

% Check four cases for prox calculation

% Case 1
if (eta+y_norm2/(2*gamma)<=0) && (y_norm <= gamma*rho)
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif (eta <= -gamma*rho^2/2) && (y_norm>gamma*rho)
    
    etaProx = 0;
    yProx = y-gamma*rho*y./y_norm;
    caseInds = 2;
    
    % Case 3
elseif eta > -gamma*rho^2/2 && y_norm > rho*eta + gamma*rho*(1+rho^2/2);% + gamma*rho*(1+rho/2)
    
    etaProx = eta + gamma*rho^2/2;
    yProx = y  - gamma*rho*y./y_norm;
    
    caseInds = 3;
    
    % Case 4
elseif eta > -gamma*rho^2/2 && y_norm <= rho*eta + gamma*rho*(1+rho^2/2);
    
    alpha = 2;
    
    if 4*gamma*eta + 2*y_norm2 > 0
        
        % Compute mu
        mu = (4/alpha^2)*sqrt(y_norm2./(gamma^2) + (4/(27*alpha^2)*((alpha*eta)/gamma + 2)^3));
        
        p_s1 = gamma + (alpha*eta)/2;
        
        rt3term = 4*y_norm/(alpha^2*gamma);
        
        if imag(mu)==0
            p_denom = p_s1 + alpha^2*gamma/8*(sign(rt3term + mu)*abs(rt3term + mu)^(1/3) ...
                + sign(rt3term - mu)*abs(rt3term - mu)^(1/3)).^2;
        else
            p_denom = p_s1 + alpha^2*gamma/8*((rt3term + mu)^(1/3) ...
                + (rt3term - mu)^(1/3)).^2;
            warning('Imaginary part is not zero')
        end
        
        p = y/p_denom;
        p_norm = sqrt(sum(p.^2));
        
        % Test polynomial equation
        %t = y_norm/p_denom
        %res = t^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*t - 8*y_norm/(alpha^2*gamma)
        
        etaProx = eta + 1/2*gamma*(alpha*p_norm^2/2);
        yProx = y-gamma*p;
        %[etaProx,yProx] = proxgTREX(eta,y,0,0,2,gamma);
        
    else
        
        etaProx = 0;
        yProx = zeros(n,1);
        
    end
    %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
    %res2 = p^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*p - 8*y_norm/(alpha^2*gamma)
    caseInds = 4;
    
else
    warning('Case not covered')
    
    etaProx = eta;
    yProx = y;
    caseInds = 5;
    
end

eta_yProx = [0;Yvec]+[etaProx;yProx];

end


%%
function [eta_yProx] = proxgsHuber(eta_y,rho,gamma,Yvec)
% Perspective of Huber in 2 dimensions

% Number of samples
n = length(Yvec);

etaVec = eta_y(1:n);

% Shift by data vector
y = (eta_y(n+1:end)-Yvec);

% Number of samples
n = length(y);

etaProx = zeros(n,1);
yProx = zeros(n,1);

% Iterate over all data points
for i=1:n
    
    % ith component
    y_i = y(i);
    
    eta_i = etaVec(i);
    
    % Check five cases for prox calculation
    yabs = abs(y_i);
    
    if (eta_i+yabs^2/(2*gamma)<=0) && (yabs <= gamma*rho)
        etaProx(i)=0;
        yProx(i)=0;
        
    elseif (eta_i <= -gamma*rho^2/2) && (yabs>gamma*rho)
        
        etaProx(i) = 0;
        yProx(i) = y_i-gamma*rho*sign(y_i);
        
    elseif eta_i > -gamma*rho^2/2 && y_i > rho*eta_i + gamma*rho*(1+rho^2/2);% + gamma*rho*(1+rho/2)
        
        etaProx(i) = eta_i + gamma*rho^2/2;
        yProx(i) = y_i - gamma*rho;
        
    elseif eta_i > -gamma*rho^2/2 && y_i < -rho*eta_i - gamma*rho*(1+rho^2/2);% - gamma*rho*(1+rho/2)
        
        etaProx(i) = eta_i + gamma*rho^2/2;
        yProx(i) = y_i + gamma*rho;
        
    elseif eta_i > -gamma*rho^2/2 && yabs <= rho*eta_i + gamma*rho*(1+rho^2/2);
        
        yNorm2 = sum(y_i.^2);
        yNorm = sqrt(yNorm2);
        alpha = 2;
        
        if 4*gamma*eta_i + 2*yNorm2 > 0
            
            % Compute mu
            mu = (4/alpha^2)*sqrt(yNorm2./(gamma^2) + (4/(27*alpha^2)*((alpha*eta_i)/gamma + 2)^3));
            
            p_s1 = gamma + (alpha*eta_i)/2;
            
            rt3term = 4*yNorm/(alpha^2*gamma);
            
            p_denom = p_s1 + alpha^2*gamma/8*(sign(rt3term + mu)*abs(rt3term + mu)^(1/3) ...
                + sign(rt3term - mu)*abs(rt3term - mu)^(1/3)).^2;
            p = y_i/p_denom;
            p_norm = sqrt(sum(p.^2));
            
            % Test polynomial equation
            %t = yNorm/p_denom
            %res = t^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*t - 8*yNorm/(alpha^2*gamma)
            
            etaProx(i) = eta_i + 1/2*gamma*(alpha*p_norm^2/2);
            yProx(i) = y_i-gamma*p;
            %[etaVec(i),yProx] = proxgTREX(eta,y,0,0,2,gamma);
            
        else
            warning('set to zeros')
            etaProx(i) = 0;
            yProx(i) = 0;
            
        end
        %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
        %res2 = p^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*p - 8*yNorm/(alpha^2*gamma)
        
    else
        warning('Case not covered')
        
        etaProx(i) = eta_i;
        yProx(i) = y_i;
        
    end
    
    
end

% % Shift by data vector and return (n + n) dimensional vector
eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];

end

%%
function eta_yProx = proxgpVapnik(eta_y,epsV,gamma,Yvec)
% Perspective of Vapnik in (n+1) dimensions

eta = eta_y(1);
y = eta_y(2:end)-Yvec;

% Dimension of the problem
n = length(y);

% Norm of y
y_norm2 = sum(y.^2);
y_norm = sqrt(y_norm2);

% Check five cases for prox calculation

% Case 1
if (eta+epsV*y_norm<=0) && (y_norm <= gamma)
    etaProx=0;
    yProx=zeros(n,1);
    % Case 2
elseif (eta<=-gamma*epsV) && (y_norm > gamma)
    
    etaProx = 0;
    yProx = y-gamma*y./y_norm;
    % Case 3
elseif eta > -gamma*epsV && y_norm > epsV*eta + gamma*(1+epsV^2);
    
    etaProx = eta + gamma*epsV;
    yProx = y - gamma*y./y_norm;
    
    % Case 4
elseif y_norm>-eta/epsV && epsV*eta<= y_norm  && y_norm <=epsV*eta + gamma*(1+epsV^2)
    
    etaProx = (eta+epsV*y_norm)/(1+epsV^2);
    yProx = epsV*(eta+epsV*y_norm)*y/(y_norm*(1+epsV^2));
    
    % Case 5
elseif eta >= 0 && y_norm <= epsV*eta;
    
    etaProx = eta;
    yProx = y;
    
else
    
    warning('Case not covered')
    
    etaProx = eta;
    yProx = y;
    
end

eta_yProx = [0;Yvec] + [etaProx;yProx];

end


%%
function [eta_yProx] = proxgsVapnik(eta_y,epsV,gamma,Yvec)
% Perspective of Vapnik in 2 dimensions

% Number of samples
n = length(Yvec);

etaVec = eta_y(1:n);

% Shift by data vector
y = (eta_y(n+1:end)-Yvec);

% Number of samples
n = length(y);

etaProx = zeros(n,1);
yProx = zeros(n,1);

% Iterate over all data points
for i=1:n
    
    % ith component
    y_i = y(i);
    
    eta_i = etaVec(i);
    
    % Check five cases for prox calculation
    yabs = abs(y_i);
    
    if (eta_i+epsV*yabs<=0) && (yabs <= gamma)
        etaProx(i)=0;
        yProx(i) = 0;
        
    elseif (eta_i<=-gamma*epsV) && (yabs > gamma)
        
        etaProx(i) = 0;
        yProx(i) = y_i-gamma*sign(y_i);
        
    elseif eta_i > -gamma*epsV && yabs > epsV*eta_i + gamma*(1+epsV^2);
        
        etaProx(i) = eta_i + gamma*epsV;
        yProx(i) = y_i - gamma*sign(y_i);
        
    elseif yabs>-eta_i/epsV && epsV*eta_i<= yabs  && yabs <= (epsV*eta_i + gamma*(1+epsV^2))
        
        etaProx(i) = (eta_i+epsV*yabs)/(1+epsV^2);
        yProx(i) = epsV*(eta_i+epsV*yabs)*y_i/(yabs*(1+epsV^2));
        
        % Case 5
    elseif eta_i >= 0 && yabs <= epsV*eta_i;
        
        etaProx(i) = eta_i;
        yProx(i) = y_i;
        
    else
        
        warning('Case not covered')
        
        etaProx(i) = eta_i;
        yProx(i) = y_i;
        
    end
end

% % Shift by data vector and return (n + n) dimensional vector
eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];

end

function eta_yProx = proxgqsBerHu(eta_y,rho,alpha,kappa,qPower,gamma,Yvec)

% Number of samples
n = length(Yvec);

etaVec = eta_y(1:n);

% Shift by data vector
y = (eta_y(n+1:end)-Yvec);


etaProx = zeros(n,1);
yProx = zeros(n,1);

for i=1:n
    
    y_i = y(i);
    
    eta_i = etaVec(i);
    
    [etaProx_i,yProx_i] = proxgqBerHu(eta_i,y_i,rho,alpha,kappa,qPower,gamma);
    etaProx(i) = etaProx_i;
    yProx(i) = yProx_i;
    
end

% % Shift by data vector and return (n + n) dimensional vector
eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];
end

% Perspective of generalized BerHu in 2D

function eta_yProx = proxgsBerHu(eta_y,rho,gamma,Yvec)
% Perspective of BerHu in 2D

% Number of samples
n = length(Yvec);

etaVec = eta_y(1:n);

% Shift by data vector
y = (eta_y(n+1:end)-Yvec);


etaProx = zeros(n,1);
yProx = zeros(n,1);

for i=1:n
    
    y_i = y(i);
    
    eta_i = etaVec(i);
    
    % Check four cases for prox calculation
    yabs = abs(y_i);
    
    %%
    if eta_i > 0 && yabs < gamma
        
        etaProx(i)=eta_i;
        yProx(i)=0;
        
        %%
    elseif yabs < gamma + rho*eta_i && yabs > gamma % eta>0
        
        etaProx(i) = eta_i;
        yProx(i) = y_i-gamma*sign(y_i);
        
        %%
    elseif 2*gamma*eta_i + rho*(yabs^2-gamma^2)>0 && yabs > gamma+eta_i*rho
        
        
        fac3 = 1;
        fac2 = 0;
        fac1 = 2/(gamma*rho^2)*(gamma+rho*eta_i)-1;
        fac0 = -2/(gamma*rho^2)*yabs;
        rootfacs = [fac3 fac2 fac1 fac0];
        cubroots = roots(rootfacs);
        
        [~,ind] = max(real(cubroots));
        
        
        %ind2 = intersect(find(imag(cubroots)==0),find(real(cubroots)>=0));
        %if ind2~=ind
        %    error('indices not the same')
        %end
        
        % Plot cubic root
        %cubroots
        
        t = cubroots(ind);
        
        p_denom = gamma+rho*(eta_i+(gamma*rho/2)*(t^2-1));
        
        pp = 1/p_denom * y_i;
        
        if abs(pp)<1
            error('Wrong p');
        end
        
        % Test polynomial equation
        
        etaProx(i) = eta_i+gamma*((rho/2)*(pp^2-1));
        yProx(i) = y_i-gamma*pp;
        
    elseif 2*gamma*eta_i + rho*max(yabs^2-gamma^2,0) <= 0
        
        etaProx(i) = 0;
        yProx(i) = 0;
        
        
    else
        warning('Case not covered')
        etaProx(i) = eta_i
        yProx(i) = y_i
    end
    
end

% % Shift by data vector and return (n + n) dimensional vector
eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];
end

% Prox for Lq objective
function etaProx_yProx = proxgLq(eta_y,alpha,qPower,gamma,Yvec)

eta = eta_y(1);
y = eta_y(2:end);

% Shifted eta
etaPrime = eta;

% Dual norm
q_s = qPower/(qPower-1);

if (abs(q_s-round(q_s)))>1e-6 && qPower~=3
    error('Dual norm of q must be integer!')
end

% Factor
delta = (alpha*(1-(1/q_s)))^(q_s-1);

% Check whether prox calculation is needed

% First summand of prox check
q_t = q_s*gamma^(q_s-1)*etaPrime;

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);
yDiffNorm = sqrt(yDiffNorm2);

testVal = (q_t + delta*(yDiffNorm.^q_s));

% Check whether prox calculation is needed

if testVal>0
    
    % Solve the q* root using matlab
    if qPower==3
        % Highest root + 1
        rootfacs = zeros(1,5);
        rootfacs(5) = 1;
        rootfacs(3) = q_s.*etaPrime/(gamma*delta);
        rootfacs(2) = q_s/(delta^2); % for q=q_s=2 we already have a component
        rootfacs(1) = -q_s*yDiffNorm./(gamma*delta^2);
    else
        q_s = round(q_s);
        
        % Highest root + 1
        rootfacs = zeros(1,2*q_s);
        rootfacs(2*q_s) = 1;
        rootfacs(q_s) = q_s.*etaPrime/(gamma*delta);
        rootfacs(2) = rootfacs(2)+q_s/(delta^2); % for q=q_s=2 we already have a component
        rootfacs(1) = -q_s*yDiffNorm./(gamma*delta^2);
 
    end
    
    % Reverse root factors (factor for highest degree first,...)
    rootfacs = rootfacs(end:-1:1);
    qroots = roots(rootfacs);
    
    % Only the real root
    t = qroots(imag(qroots)==0);
    
    if sum(t>=0)>1
        error('Too many positive real roots')
    end
    
    % Take the largest real root
    t = max(t);
    
    if qPower==3
        t = sqrt(t);
    end
    
    % p_denom  = gamma+delta*(eta+gamma*delta*(t.^q_s)/q_s)*t.^(q_s-2);
    % p = y/(p_denom);
    
    % Simplified prox computation
    p = (yDiff/yDiffNorm) * t;
    
    etaProx_yProx = [eta + gamma*delta*(t.^q_s)/q_s;y-gamma*p];
    
else
    
    etaProx_yProx = [0;Yvec];
    
end
end

% Prox for element-wise Lq objective
function eta_yProx = proxgsLq(eta_y,alpha,qPower,gamma,Yvec)

% Number of samples
n = length(Yvec);

etaVec = eta_y(1:n);

% Shift by data vector
y = (eta_y(n+1:end)-Yvec);

etaProx = zeros(n,1);
yProx = zeros(n,1);

% Dual norm
q_s = qPower/(qPower-1);

if (abs(q_s-round(q_s)))>1e-6 && qPower~=3
    error('Dual norm of q must be integer!')
end

% Factor
delta = (alpha*(1-(1/q_s)))^(q_s-1);

for i=1:n
    
    y_i = y(i);
    eta_i = etaVec(i);
    
    % First summand of prox check
    q_t = q_s*gamma^(q_s-1)*eta_i;
    
    % Check cases for prox calculation
    yabs = abs(y_i);
    
    testVal = (q_t + delta*(yabs.^q_s));
    
    
    %disp(['Delta value: ',num2str(delta)])
    %disp(['qt value: ',num2str(q_t)])
    %disp(['Test value: ',num2str(testVal)])
    
    % Check whether prox calculation is needed
    if testVal>0
        
        % Solve the q* root using matlab
        if qPower==3
            % Highest root + 1
            rootfacs = zeros(1,5);
            rootfacs(5) = 1;
            rootfacs(3) = q_s.*eta_i/(gamma*delta);
            rootfacs(2) = q_s/(delta^2); % for q=q_s=2 we already have a component
            rootfacs(1) = -q_s*yabs./(gamma*delta^2);
        else
            q_s = round(q_s);
            
            % Highest root + 1
            rootfacs = zeros(1,2*q_s);
            rootfacs(2*q_s) = 1;
            rootfacs(q_s) = q_s.*eta_i/(gamma*delta);
            rootfacs(2) = rootfacs(2)+q_s/(delta^2); % for q=q_s=2 we already have a component
            rootfacs(1) = -q_s*yabs./(gamma*delta^2);
            
        end
        
        % Reverse root factors (factor for highest degree first,...)
        rootfacs = rootfacs(end:-1:1);
        qroots = roots(rootfacs);
        
        % Only the real root
        t = qroots(imag(qroots)==0);
        
        if sum(t>=0)>1
            error('Too many positive real roots')
        end
        
        % Take the largest real root
        t = max(t);
        
        if qPower==3
            t = sqrt(t);
        end
        
        % p_denom  = gamma+delta*(eta+gamma*delta*(t.^q_s)/q_s)*t.^(q_s-2);
        % p = y/(p_denom);
        
        % Simplified prox computation
        p = (y_i/yabs) * t;
        
        etaProx(i) = eta_i + gamma*delta*(t.^q_s)/q_s;
        yProx(i) = y(i)-gamma*p;
        
    else
        
        etaProx(i) = eta_i;
        yProx(i) = y(i);
        
    end
end

eta_yProx = [zeros(n,1);Yvec]+[etaProx;yProx];
end

%%% Concomitant Lasso
function etaProx_yProx = proxgpL2(eta_y,alpha,gamma,Yvec)

eta = eta_y(1);
y = eta_y(2:end);

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);

% Check whether prox calculation is needed
if (4*gamma*eta+alpha*yDiffNorm2)>0
    
    %disp(['Hello 2 with eta=',num2str(eta)])
    
    yDiffNorm = sqrt(yDiffNorm2);
    
    % Compute mu
    mu = (4/alpha^2)*sqrt(yDiffNorm2/(gamma.^2) + 32/(27*alpha^2)*((alpha*eta)/(2*gamma)+1)^3);
    
    % Compute root explicitly
    p_s1 = gamma + (alpha*eta)/2;
    
    rt3term = 4*yDiffNorm/(alpha^2*gamma);
    
    p_denom = p_s1 + ((alpha^2*gamma)/8)*(sign(rt3term + mu)*(abs(rt3term + mu))^(1/3) ...
        + sign(rt3term - mu)*(abs(rt3term - mu))^(1/3)).^2;
    %yDiffNorm
    p = yDiff./(p_denom);
    
    % Root checking
    % t = yDiffNorm/p_denom;
    % t^3+4*(alpha*eta+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
    
    % Compute polynomial
    etaProx_yProx = [eta + alpha*gamma*sum((p).^2)/4;y-gamma*p];
    
else
    
    etaProx = 0;
    yProx = Yvec;
    etaProx_yProx = [etaProx;yProx];
    
end
end



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
