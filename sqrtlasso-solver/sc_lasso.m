function [betaSCMat, sigmaSCMat,funSCMat,out] = sc_lasso(Xmat, Yvec, inopts)
% Smooth Concomittant Lasso (Ndiaye et al, 2016)
% Approximate solver using coordinate descent (soft thresholding)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: SC-LASSO solutions: betaSCMat,sigmaSCMat,funSCMat

[n, p] = size(Xmat);

% Options for prox
defopts.abstol = 1e-4;          % Primal-dual gap tolerance
defopts.gapit = 1e1;            % Interval for checking primal dual gap
defopts.maxit = p*2e2;          % p*1e2 Number of iteration in Coordinate descent
defopts.sigma_0_fac = 1e-1;     % Data-dependent scaling of lower bound on sigma
defopts.lenPath = 200;          % 200 Length of lambda path
defopts.delta = 2;              % 2 Constant for geometric sequence spacing of the regularization path
defopts.activeSet = 1:p;        % Active set of variables
defopts.lamVec = [];            % User-defined lambda values

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end


% Active set on which SC-LASSO is computed in the beginning (will be
% modified by the screening rule)
actInds = opts.activeSet;
nA = length(actInds);


if ~isempty(opts.lamVec)
    opts.lenPath = length(opts.lamVec);
    lamVec = opts.lamVec;
    
    % Length of the regularization path
    nLams = opts.lenPath; 
else
    
    
    % Length of the regularization path
    nLams = opts.lenPath;
    
    % Data-dependent lambda upper bound
    lam_max = max(abs(Xmat'*Yvec))./(sqrt(sum(Yvec.^2))*sqrt(n-1));
    
    % Initial lambda
    lam_0 = lam_max;
    
    % Geometric sequence for lambda path
    if nLams>1
        delta = opts.delta;
        lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1))
    end
end

% Initial beta
beta_0 = zeros(p,1);

% Create storage for 2 nA beta vectors, function values and timings
betaSCMat = zeros(p,nLams);
sigmaSCMat = zeros(1,nLams);
funSCMat = zeros(1, nLams);
dualityCell = cell(nLams,1);
runTimes = zeros(1, nLams);

% Algorithm stopping criteria
maxit = opts.maxit;
abstol = opts.abstol;

% Interval for gap computation
gapit = opts.gapit;

% Initial sigma
sigma_i = sqrt(sum(Yvec.^2))/sqrt(n);

% Lower bound on sigma for SC-Lasso
sigma_0 = opts.sigma_0_fac * sigma_i;

% Starting vector
beta_i = beta_0;

% Precompute properties of the data ONCE
XY = Xmat'*Yvec;
XX = Xmat'*Xmat;
sum_X2 = sum(Xmat.^2);

% Loop over lambda path
for i = 1:nLams
    
    gapVec = [];
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i)
    
    % Variables needed
    nlam_i = n*lam_i;
    sqrt_n = sqrt(n);
    
    % Iterations for coordinate descent
    for k=1:maxit
        
        if ~mod(k,gapit)
            
            % Primal-dual gap computation
            primSol = 1/(2*n*sigma_i)*sum((Yvec-Xmat*beta_i).^2) + sigma_i/2 + lam_i*sum(abs(beta_i));
            % Equation (9) in paper
            dualDenom1 = (n*lam_i*sigma_0);
            dualDenom2 = max(Xmat'*(Yvec-Xmat*beta_i));
            dualDenom3 = lam_i*sqrt(n)*sqrt(sum((Yvec-Xmat*beta_i).^2));
            dualDenom = max([dualDenom1,dualDenom2,dualDenom3]);
            
            dualTheta = (Yvec-Xmat*beta_i)./dualDenom;
            
            dualSol = Yvec'*lam_i*dualTheta + sigma_0*(1/2-lam_i^2*n/2*sum(dualTheta.^2));
            
            % Primal dual gap
            G_sig_lam = (primSol-dualSol);
            
            % Store duality gap
            gapVec = [gapVec,G_sig_lam];
            
            if G_sig_lam<abstol
                %disp('Primal-dual gap reached')
                break
            end
            
            % Screening rule
            r_safe = sqrt((2*G_sig_lam)./(lam_i^2*sigma_0*n));
            safeSphereVals = (abs(Xmat'*dualTheta) + r_safe*sqrt(sum_X2'));
            actInds = find(safeSphereVals>1)';
            if length(actInds)<p
                disp(['Screening active: p_act=',num2str(length(actInds))])
            end
        end
        
        % Initialize matrix vector muliply
        Xmat_beta = Xmat*diag(beta_i);
        
        % Coordinate descent over all active variables
        for j=actInds
            sum_j = sum_X2(j);
            proxArg = beta_i(j)-(XX(j,:)*beta_i-XY(j))/sum_j;
            gamma = sigma_i*nlam_i/sum_j;
            beta_i(j) = proxST(proxArg,gamma);
            Xmat_beta(:,j) = Xmat(:,j)*beta_i(j);
            % Scaled lasso style
            sigma_i = max(sigma_0,sqrt(sum((Yvec-sum(Xmat_beta,2)).^2))/sqrt_n);
            %% TREX style
            %%sigma_i = max(sigma_0,max((XY-XX*beta_i)/(lam_i*n))); % TREX style
            
        end
        %figure(23)
        %plot(k,sigma_i,'b.','MarkerSize',50)
        %hold on
        %grid on
        %drawnow
    end
    
    betaSCMat(:,i) = beta_i;
    sigmaSCMat(1,i) = sigma_i;
    funSCMat(1,i) = 1/(2*n*sigma_i)*sum((Yvec-Xmat*beta_i).^2) + sigma_i/2 + lam_i*sum(abs(beta_i));
    runTimes(1,i) = toc;
    dualityCell{i} = gapVec;
end

out.opts = opts;
out.runTimes = runTimes;
out.dualityCell = dualityCell;
out.X = Xmat;
out.Y = Yvec;
out.lamPath = lamVec;
out.sigma_0 = sigma_0;
end

% Soft thresholing
function betaT = proxST(b,gamma)
betaT = sign(b).*max(abs(b) - gamma,0);
end

function [funVal,L1beta] = objTREX_SA(betaVec,X,Y,c)
% TREX standalone objective function
% funVal: function value at current betaVec
% L1beta: L1 norm of betaVec

% Global data vector X, response Y, normalizing constant normConst, and
% q-norm approximation to inf-norm

[n,p]= size(X);

funVal = sum((Y-X*betaVec).^2)/(c*max(abs(X'*(Y-X*betaVec))));

L1beta = sum(abs(betaVec));

funVal = funVal+L1beta;
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
