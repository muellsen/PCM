function betaT = proxL1(b,lam,gam,weights)

% Proximity operator of L1 norm

% Regularization parameter 
lam_gam = lam*gam;

% Soft thresholding of the predictors
betaT = sign(b).*max(abs(b) - lam_gam.*weights,0);
