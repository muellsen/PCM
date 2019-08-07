%% input:
% x: n by p design matrix
% y: response of length n
% constr: linear constraint, for compositional data, constr=ones(p,1)/sqrt(p); for no constraint, constr=zeros(p,1)
% level: level of confidence interval
% tol: tolerance level in the computation
% maxiter: the maximum number of iteration in solving the regularized estimator using coordinate descent method

x = normrnd(0,1,100,120);
bet = [1,0.8,0.6,-0.5,-0.7,-1.2,zeros(1,120-6)]';
constr = ones(120,1);
y = x*bet + normrnd(0,0.5);
level = 0.95;
tol = 1e-8;
maxiter = 200;



%% compute the projection matrix according to the constraints
[constr2,S,V] = svd(constr);
constr = constr2(:,1:size(S,2));
Pc = constr*constr';
x = x-x*Pc;

%% get regularized estimator
method = 1; % using coordinate descent
method = 2; % using CVX
[res.bet_n res.int res.sigma res.lam0]=sslr(y, x, constr, tol, maxiter, method);

%% get confidence interval
% for faster computation, we suggest you use CVX when p>2*n
[res.bet_u res.CI_l res.CI_u res.M]=cdmm_ci(y, x, constr, res.bet_n, res.int, res.sigma, res.lam0, level, tol, maxiter); % using coodinate descent
[res.bet_u res.CI_l res.CI_u res.M]=cvx_ci(y, x, constr, res.bet_n, res.int, res.sigma, res.lam0, level, tol); % using CVX

%% output
% bet_n: regularized estimator
% int: intercept
% sigma: estimated noise level
% bet_u: debiased estimator
% [CI_l, CI_u]: confidence interval

