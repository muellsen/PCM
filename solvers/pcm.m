function [betaPCMMat, sigmaPCMMat,funPCMMat,out] = pcm(Xmat, Yvec, inopts)
% Penalized Concomitant M-estimators (Combettes et al, 2017)
% Generic proximal solver for penalized concomitant estimators
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: PCM solutions: betaSCMat,sigmaSCMat,funSCMat

[n, p] = size(Xmat);

% Options for regularization path
defopts.lenPath = 200;          % Length of lambda path
defopts.delta = 2;              % Constant for geometric sequence spacing of the regularization path
defopts.activeSet = 1:p;        % Active set of variables
defopts.lambda = [];            % Single lambda value
defopts.plotting = 0;           % Plotting of trace
defopts.verbose = 0;

% Options for proximal scheme
defopts.abstol = 1e-12;         % DR tolerance
defopts.maxit = p*1e2;          % p*1e2 Number of iteration in DR
defopts.rho = 1.345;            % rho value for Huber function
defopts.epsV = 0.1;             % eps value for Vapnik function    
defopts.gamma = 20;             % gamma in DR

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

% Active set on which SC-LASSO is computed in the beginning (will be
% modified by the screening rule)
actInds = opts.activeSet;
nA = length(actInds);

% Length of the regularization path
nLams = opts.lenPath;

% Create storage for 2 nA beta vectors, function values and timings
betaPCMMat = zeros(p,nLams);
sigmaPCMMat = zeros(1,nLams);
funPCMMat = zeros(1, nLams);
runTimes = zeros(1, nLams);
betaOrbitCell = cell(1,nLams);

% Algorithm stopping criteria
maxit = opts.maxit;
abstol = opts.abstol;

% DR gamma
gamma = opts.gamma;

% Huber rho
rho = opts.rho

% Vapnik eps
epsV = opts.epsV;

% Data-dependent lambda upper bound
%lam_max = max(Xmat'*Yvec)./(sqrt(sum(Yvec.^2))*sqrt(n));
%lam_max = max(Xmat'*Yvec)./(sqrt(sum(Yvec.^2)));

% For Huber penalty 
rhoInds = (abs(Yvec)>rho);
nrhoInds = (abs(Yvec)<=rho);
lam_max = max(Xmat'*Yvec)./sum([abs(Yvec(rhoInds));(Yvec(nrhoInds)).^2/2])

% Initial beta (solution vector)
beta0 = zeros(p,1);

% Initial lambda
lam_0 = lam_max;

% Geometric sequence for lambda path if not set externally
if nLams>1 && isempty(lamVec)
    delta = opts.delta;
    lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1));
end

%% Reverse trace
%lamVec = lamVec(end:-1:1);



% Loop over regularization path
for i = nLams:-1:1
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i);
    
    % Douglas-Rachford for Huber
    [beta_K,betaOrbit,M,R] = dougRach(Yvec,Xmat,beta0,sigma0,epsV,rho,gamma,lam_i,opts);
    
    betaPCMMat(:,i) = beta_K(2:end,1);
    sigmaPCMMat(1,i) = beta_K(1);
    betaOrbitCell{1,i} = betaOrbit;
    runTimes(1,i) = toc;
    
    % Exit when Vapnik is violated
%     if norm(beta_K(2:end,1))<1e-1
%         disp(['Upper bound reached at $\lambda$: ',num2str(lam_i)])
%         break
%     end
    
    %Warm start
    sigma0 = beta_K(1);
    beta0 = beta_K(2:end,1);
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
function [betaPCM,betaOrbit,M,R] = dougRach(Yvec,Xmat,beta0,sigma0,epsV,rho,gamma,lam,dropts)

maxit = dropts.maxit;
abstol = dropts.abstol;
mu_k = dropts.dr_mu;

[n,p] = size(Xmat);

% n+1 x p matrix
M = zeros(n+1,p+1);
M(1,1) = 1;
%M(2:n+1,1) = Yvec;
M(2:n+1,2:p+1) = Xmat;

% Matrices
Id_n = eye(n+1);
R = M'/(Id_n+M*M');

% Initialize DR sequences

y0 = Yvec;

x_k = [sigma0;beta0]; % p+1 vector
y_k = [sigma0;y0];    % n+1 vector

% Help variable
z_k1 = [sigma0;beta0];

betaOrbit = zeros(p+1,maxit);

% Main DR loop
for kk=1:maxit
    
    q_k = M*x_k - y_k;
    w_k = x_k - R*q_k;
    r_k = M*w_k;
    
    % Soft thresholding
    z_k = proxST(2*w_k-x_k,lam*gamma);

    % Ridge regularization
    %z_k = proxL2(2*w_k-x_k,lam*gamma);
    
    % Data fitting
    t_k = proxgsHuber(2*r_k-y_k,rho,gamma,Yvec);
    
    %t_k = proxgpHuber(2*r_k-y_k,rho,gamma,Yvec);
    %t_k = proxgVapnik(2*r_k-y_k,epsV,gamma,Yvec);

    %t_k = proxL2(2*c_k-y_k,Yvec,gamma);
    %t_k = 2*c_k-y_k;
    
    % currSigma = max(0.1,1/sqrt(n)*sqrt(sum((Xmat*b_k(2:end)-Yvec).^2)))
    
    x_k1 = x_k + mu_k*(z_k-w_k);
    y_k1 = y_k + mu_k*(t_k-r_k);
    
    %norm(x_k1-x_k)
    %norm(z_k1-z_k)
    %abs(norm(x_k1-x_k)-norm(z_k1-z_k))
    
    % x_k-z_k
    % y_k-t_k
    xx_k = [w_k;r_k];
    zz_k = [z_k;t_k];
    
    norm_xz = norm(xx_k-zz_k);
    
    if ~mod(kk,2e3)
        disp(['Error at iteration: ',num2str(kk),': ',num2str(norm_xz)]);
    end
        
    if (norm_xz<abstol) && (kk>1e2)
        x_k = x_k1;
        y_k = y_k1;
        
        disp(['Early convergence at iteration: ',num2str(kk)])
        break
    end
    
    z_k1 = z_k;
    
    betaOrbit(:,kk) = w_k;
    x_k = x_k1;
    y_k = y_k1;
    
    % Adapting mu
    %mu_k = (1-1e-3)*mu_k;
    
    if dropts.plotting
        
        if ~mod(kk,2e1)
            semilogy([K_old,kk],abs([w_old,w_k]),'k-','LineWidth',5)
            grid on
            hold on
            drawnow
            w_old = w_k;
            K_old = kk;
        end
    end
end

% Final estimate (sigma is first component, beta follows)
betaPCM = w_k;

end


function betaT = proxST(b,lam_gam)
% Soft thresholding operation

betaT = b;

% Do not threshold first entry (contains sigma)
betaT(2:end) = sign(b(2:end)).*max(abs(b(2:end)) - lam_gam,0);

end

function eta_yProx = proxL2(eta_y,gamma)
% Proximity operator for L2

etaProx = eta_y(1);
y = eta_y(2:end);

yProx = (1-gamma./norm(y))*y;

eta_yProx = [etaProx;yProx];

end

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



function [eta_yProx] = proxgsHuber(eta_y,rho,gamma,Yvec)
% Perspective of Huber in 2 dimensions

eta = eta_y(1);
y = eta_y(2:end)-Yvec;

% Number of samples 
n = length(y);

etaVec = zeros(n,1);
yProx = zeros(n,1);

% Iterate over all data points
for i=1:n
    
    % ith component
    y_i = y(i);
    
    % Check five cases for prox calculation
    yabs = abs(y_i);
        
    if (eta+yabs^2/(2*gamma)<=0) && (yabs <= gamma*rho)
        etaVec(i)=0;
        yProx(i)=0;      
        
    elseif (eta <= -gamma*rho^2/2) && (yabs>gamma*rho)
        
        etaVec(i) = 0;
        yProx(i) = y_i-gamma*rho*sign(y_i);
        
    elseif eta > -gamma*rho^2/2 && y_i > rho*eta + gamma*rho*(1+rho^2/2);% + gamma*rho*(1+rho/2)
        
        etaVec(i) = eta + gamma*rho^2/2;
        yProx(i) = y_i - gamma*rho;
        
    elseif eta > -gamma*rho^2/2 && y_i < -rho*eta - gamma*rho*(1+rho^2/2);% - gamma*rho*(1+rho/2)
        
        etaVec(i) = eta + gamma*rho^2/2;
        yProx(i) = y_i + gamma*rho;
                
    elseif eta > -gamma*rho^2/2 && yabs <= rho*eta + gamma*rho*(1+rho^2/2);
                
        yNorm2 = sum(y_i.^2);
        yNorm = sqrt(yNorm2);
        alpha = 2;
        
        if 4*gamma*eta + 2*yNorm2 > 0
            
            % Compute mu
            mu = (4/alpha^2)*sqrt(yNorm2./(gamma^2) + (4/(27*alpha^2)*((alpha*eta)/gamma + 2)^3));
            
            p_s1 = gamma + (alpha*eta)/2;
            
            rt3term = 4*yNorm/(alpha^2*gamma);
            
            p_denom = p_s1 + alpha^2*gamma/8*(sign(rt3term + mu)*abs(rt3term + mu)^(1/3) ...
                + sign(rt3term - mu)*abs(rt3term - mu)^(1/3)).^2;
            p = y_i/p_denom;
            p_norm = sqrt(sum(p.^2));
            
            % Test polynomial equation
            %t = yNorm/p_denom
            %res = t^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*t - 8*yNorm/(alpha^2*gamma)
            
            etaVec(i) = eta + 1/2*gamma*(alpha*p_norm^2/2);
            yProx(i) = y_i-gamma*p;
            %[etaVec(i),yProx] = proxgTREX(eta,y,0,0,2,gamma);
            
        else
            
            etaVec(i) = 0;
            yProx(i) = 0;
            
        end
        %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
        %res2 = p^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*p - 8*yNorm/(alpha^2*gamma)
            
    else
        warning('Case not covered')
        
        etaVec(i) = eta;
        yProx(i) = y_i;
        
    end
    
    
end

% Projection onto the mean of the concomitant
etaProx = mean(etaVec);

eta_yProx = [0;Yvec]+[etaProx;yProx];

end


function eta_yProx = proxgVapnik(eta_y,epsV,gamma,Yvec)
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
