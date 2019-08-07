function [betaVec, sigmaVec, numIter]=sslr_lam(y, x, constr, lam0, tol,varargin)
% Iterative solver for concomitant estimation of sigma and beta
% for fixed lambda

% Varargin may contain a beta vector for warm starts

[n p] = size(x);

% Initial sigma values
sigma = 1;
sigma_2 = 1;
sigma_s = 0.5;

% Maximum iteration for solver
maxIter = 1e3;

% Numerical tolerance for zeros
numZTol = 1e-5;
i=1;

while (abs(sigma-sigma_s)>tol & i<maxIter)
    i=i+1;
    sigma = (sigma_s + sigma_2)/2;
    lam = sigma*lam0;
    %betaV = bycvx(y, x, constr, lam);
    betaV = lasso_constr(y, x, constr, lam, 1e-6, 1e3,varargin);
    %s = sum(abs(betaV)>numZTol);
    %s = min(s, n-1);
    s=-1;
    sigma_s = norm(y-x*betaV)/sqrt(n-s-1);
    sigma_2 = sigma;
end

sigma = (sigma_s + sigma_2)/2;

sigmaVec = sigma*ones(n,1);
betaVec = betaV;
numIter = i-1;

