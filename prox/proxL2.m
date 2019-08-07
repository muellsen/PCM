function betaT = proxL2(b,lam,gam)
% Proximity operator for L2 norm

% Regularization parameter 
lam_gam = lam*gam;

b_norm = (sqrt(sum(b.^2)));

% Proximity operator
if b_norm>lam_gam
    betaT = (1-lam_gam./b_norm)*b;
else
    betaT = 0*b;
end