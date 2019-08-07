function [betaSCMat, sigmaSCMat,funSCMat,out] = huber_cvx(Xmat, Yvec, inopts)
% Concomitant Huber Lasso (Lambert-Lacroix and Zwald, 2011)
% CVX solver
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: Huber solutions: betaSCMat,sigmaSCMat,funSCMat

[n, p] = size(Xmat);

% Options for prox
defopts.rho = 1.345;            % Huber parameter rho
defopts.lenPath = 200;          % 200 Length of lambda path
defopts.sigma_0_fac = 1e-6;     % Data-dependent scaling of lower bound on sigma
defopts.adaptLasso = 0;         % Indicator whether to use adaptive Lasso
defopts.delta = 2;              % 2 Constant for geometric sequence spacing of the regularization path
defopts.lamPath = [];            % Single lambda value or path

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

if ~isempty(opts.lamPath)
    opts.lenPath = length(opts.lamPath);
    lamVec = opts.lamPath;
else
    lamVec = [];
end

% Length of the regularization path
nLams = opts.lenPath;

% Create storage for 2 nA beta vectors, function values and timings
betaSCMat = zeros(p,nLams);
sigmaSCMat = zeros(1,nLams);
funSCMat = zeros(1, nLams);
runTimes = zeros(1, nLams);

% Huber parameter rho
rho = opts.rho;

% Data-dependent lambda upper bound
% For Huber penalty
rhoInds = (abs(Yvec)>rho);
nrhoInds = (abs(Yvec)<=rho);
lam_max = max(Xmat'*Yvec)./sum([abs(Yvec(rhoInds));(Yvec(nrhoInds)).^2/2]);

% Geometric sequence for lambda path
if nLams>1 && isempty(lamVec)
    delta = opts.delta;
    lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1));
end

% Adaptive Lasso use
adaptLasso = opts.adaptLasso;

% Initial sigma
sigma_i = sqrt(sum(Yvec.^2))/sqrt(n);

% Lower bound on concomitant for Huber
s_0 = opts.sigma_0_fac * sigma_i;

% Checking whether to use adaptive lasso or not
if adaptLasso
    
    % cvx solution (OLS first)
    cvx_begin quiet
    variables betaUNP(p) s_i v(n);
    minimize (n*s_i+quad_over_lin(Yvec-Xmat*betaUNP-v,s_i)+2*rho*norm(v,1))
    subject to
    %s_i>=s_0
    s_i>0;
    cvx_end
    
    % Take absolute value as weights
    betaUNP = abs(betaUNP);
else
    % Non-adaptive Lasso
    betaUNP=ones(p,1);
end


% Loop over lambda path
for i = 1:nLams
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i);
    
    % Scaling needed for lambda path
    %snlam_i = rho*n*lam_i;
    snlam_i = lam_i;
    
    % Moreau-Yoshida formulation of concomitant Huber
    cvx_begin quiet
    variables beta_i(p) s_i v(n);
    minimize (n*s_i+quad_over_lin(Yvec-Xmat*beta_i-v,s_i)+2*rho*norm(v,1) + snlam_i*norm(beta_i./betaUNP,1))
    subject to
    %s_i>=s_0
    s_i>0;
    cvx_end
    
    betaSCMat(:,i) = beta_i;
    sigmaSCMat(1,i) = s_i;
    funSCMat(1,i) = cvx_optval;
    runTimes(1,i) = toc;
end

out.opts = opts;
out.runTimes = runTimes;
out.X = Xmat;
out.Y = Yvec;
out.lamPath = lamVec;
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
