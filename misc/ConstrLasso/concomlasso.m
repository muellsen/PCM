function [betaConLMat, sigmaConLMat,out] = concomlasso(Xmat, Yvec, inopts)
% Log-contrast constrained scaled Lasso based on Shi et al. 2018 implementation
%
% Generic cvx solver for penalized constrained concomitant Lasso
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: PCM solutions: betaPCMMat, sigmaPCMMat,funPCMMat,out

[n, p] = size(Xmat);

% Options for simple regularization term on beta vector
defopts.Ceq = ones(1,p);            % Matrix needed for log-contrast model constraints

% Options for regularization path
defopts.lenPath = 200;              % Length of lambda path
defopts.delta = 2;                  % Constant for geometric sequence spacing of the regularization path
defopts.activeSet = 1:p;            % Active set of variables
defopts.lamPath = [];               % User input lambda values
defopts.warmstart = 1;              % Option for warmstart over the lambda path

% Options for Douglas Rachford scheme
defopts.abstol = 1e-4;              % sigma tolerance


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
    lam_max = max(abs(Xmat'*Yvec))./(sqrt(sum(Yvec.^2))*sqrt(n-1));
    
    % Geometric sequence for lambda path
    delta = opts.delta;
    lamVec = lam_max*10.^(-delta*((1:nLams)-1)./(nLams-1));
    
end

% Add path back to opts structure (used in DR for warmstart checks)
opts.lamPath = lamVec;

% Log-contrast model parameters solving Ceq*b = 0;
Ceq = opts.Ceq;            

% Tolerance on sigma convergence
tol = opts.abstol;

% Temporary use of full concomitant vector independent of group identity
betaConLMat = zeros(p,nLams);
sigmaConLMat = zeros(n,nLams);

bVecInit = zeros(p,1);

% Create storage for function values and timings
runTimes = zeros(1, nLams);

% Loop over regularization path
for i = 1:nLams
    
    tic
    % Current regularization parameter
    lam_i = lamVec(i);
    
    % Iterative heuristic solver
    [bVec,sVec] = sslr_lam(Yvec, Xmat, Ceq, lam_i, tol,bVecInit);
    
    betaConLMat(:,i) = bVec;
    sigmaConLMat(:,i) = sVec;
    
    runTimes(1,i) = toc;
    
    disp(['Convergence after: ',num2str(runTimes(1,i)),' seconds'])
    
    % Break if lower bound on lambda is reached
    if sum(abs(bVec)>1e-2)>=n
        remLen = (nLams-i);
        betaConLMat(:,i+1:end)=repmat(bVec,1,remLen);
        sigmaConLMat(:,i+1:end)=repmat(sVec,1,remLen);
        runTimes(1,i+1:end) = 0;
        disp('Lower bound on lambda reached')
        break;
    end  
    
    bVecInit = bVec;
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
