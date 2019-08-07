function [bet_u, CI_l, CI_u, cov_bet]=cdmm_ci(y, x, constr, bet_n, int, sigma, lam0, level, tol, maxiter)
% compute debiased estimator and CI using coordinate descent
tol = 0.001;
maxiter = 5000;
[n p] = size(x);

gam0 = lam0/3; % for simulation

eyes = eye(p);

loc_x = mean(x);
x2 = x - ones(n,1)*loc_x;
Sig = x2'*x2/n;
Sig2 = Sig - diag(diag(Sig));

Q = eyes - constr*constr';
M = zeros(p,p);
for i = 1:p
    disp(i)
    gam = gam0/2;
    while (gam<0.5)
        gam = gam*2;
        mi = ones(p,1);
        mi0 = zeros(p,1);
        iter = 1;
        while (sum(abs(mi-mi0))>tol & iter<maxiter)
            %disp(sum(abs(mi-mi0)))
            mi0 = mi;
            for j = 1:p
                v = -Sig2(j,:)*mi + Q(j,i);
                mi(j) = sign(v) * max(0, abs(v)-gam) / Sig(j,j);
            end
            iter = iter + 1;
        end
        if (iter<maxiter)
            break
        end
    end
    disp(gam)
    M(i,:) = mi;
end

M = Q*M;
bet_u = bet_n + M*x2'*(y-x2*bet_n-int)/n;
cov_bet = sigma^2*M*Sig*M'/n;
V = diag(cov_bet);
width = norminv(1-(1-level)/2)*sqrt(V);
CI_u = bet_u + width;
CI_l = bet_u - width;



