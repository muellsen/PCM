function betaT = proxL1slin(b,lam,gam,weights)
% Proximity operator of L1 norm with sum constraints (lp solution)

% Regularization parameter
lam_gam = lam*gam*weights;

% Initialize betaT
betaT = b;

% Find constraint and unconstrained predictors
cInds = find(lam_gam~=0);
ucInds = find(lam_gam==0);

% Call function below (with numerical zero (1e-12) because of bug)
temp1 = prox_l1sumlin( b(cInds), lam_gam(cInds(1)));

% Unconstrained L1 on uc variables
temp2 = proxL1( b(ucInds),lam,gam,ones(length(ucInds),1));

% Stitch prox together
betaT(cInds) = temp1;
betaT(ucInds) = temp2;


function b_prox = prox_l1sumlin( bVec, lam_s)

d = length(bVec);

c_lam = ones(2*d,1);

A_pos = diag(-ones(2*n,1));
b_pos = zeros(2*n,1);

A_lb = [-diag(ones(n,1)),diag(ones(n,1))];
A_ub = [diag(ones(n,1)),-diag(ones(n,1))];

A = [A_pos;A_lb;A_ub];
b = [b_pos;b_bb];

%c_lam = [ones(n,1);-ones(n,1)];

%c_lam(currInds) = 1;
%c_lam(currInds+n) = 1;

b_prox = linprog(c_lam,A,b,full(S_ex),beq,[],[],[],linopts);




