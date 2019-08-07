function [tau,bT,gam_lam] = proxgL1(tau,b,gam_lam)

% Proximity operator of L1 norm
bT = sign(b).*max(abs(b) - gam_lam,0);