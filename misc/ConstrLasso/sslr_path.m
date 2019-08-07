function [bet, int, sigma, lam0,lamVec,betMat]=sslr_path(y, x, constr, tol, tol_sig,maxiter, method)
% compute the regularized/biased estimator and noise level
[n p] = size(x);

loc_y = mean(y);
loc_x = mean(x);
sca_x = std(x);
b2 = 1./sca_x;
c = diag(b2);
x2 = (x - ones(n,1)*loc_x)*c;
constr2 = constr'*c;
y2 = y - loc_y;


options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = 0.5*sqrt(2/n)*norminv(1-kk/p);
sigma = 1;
sigma_2 = 1;
sigma_s = 0.5;
i=1;

lamVec = [];
betMat = [];

if (method==2)
    while (abs(sigma-sigma_s)>tol_sig & i<500)
        i=i+1;
        sigma = (sigma_s + sigma_2)/2;
        lam = sigma*lam0;
        lamVec = [lamVec,lam];
        bet2 = bycvx(y2, x2, constr2, lam);
        betMat = [betMat,bet2];
        s = sum(abs(bet2)>0.001);
        s = min(s, n-1);
        sigma_s = norm(y2-x2*bet2)/sqrt(n-s-1);
        sigma_2 = sigma;
    end
else
    while (abs(sigma-sigma_s)>tol_sig & i<500)
        i=i+1;
        sigma = (sigma_s + sigma_2)/2;
        lam = sigma*lam0;
        lamVec = [lamVec,lam];
        bet2 = lasso_constr(y2, x2, constr2, lam, tol, maxiter);
        betMat = [betMat,bet2];
        s = sum(abs(bet2)>0.001);
        s = min(s, n-1);
        sigma_s = norm(y2-x2*bet2)/sqrt(n-s-1);
        sigma_2 = sigma;
    end
end
if i==50 disp(0);end;

sigma = (sigma_s + sigma_2)/2;
bet = c*bet2;
int = loc_y - loc_x*bet;