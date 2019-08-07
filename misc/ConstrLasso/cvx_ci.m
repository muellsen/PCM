function [bet_u, CI_l, CI_u, M]=cvx_ci(y, x, constr, bet_n, int, sigma, lam0, level, tol)
% compute debiased estimator and CI using CVX
[n p] = size(x);

gam = lam0/3; % the value used for our simulation, chosen by user

eyes = eye(p);

Sig = x'*x/n;
[U_sig,S_sig,V_sig] = svd(Sig);
S_sig = diag(S_sig);
eps = S_sig(size(S_sig,1)-1)/1000;
Sig_eps = Sig + eps*eyes;

Q = eyes - constr*constr';
b1 = ones(p,1)*gam;
M = zeros(p,p);
for i = 1:p
    disp(i)
    b2 = Q(:,i);
    l = b2 - b1;
    u = b2 + b1;

    cvx_begin quiet
    cvx_precision low
        variable m(p)
        minimize( m'*Sig_eps*m)
        subject to
        l <= Sig*m <= u
    cvx_end    
    M(i,:) = m;

    if(~(strcmp(cvx_status,'Solved')|strcmp(cvx_status,'Inaccurate/Solved')))
        disp('first_failed')
        cvx_begin quiet
         cvx_precision low
            variables gam2 m(p)
            minimize(gam2)
            subject to
            b2 - ones(p,1)*gam2 <= Sig*m <= b2 + ones(p,1)*gam2
        cvx_end
        M(i,:) = m; 
        l2 = b2 - ones(p,1)*max(1.1*gam2,1.2*gam);
        u2 = b2 + ones(p,1)*max(1.1*gam2,1.2*gam);

        cvx_begin quiet
        cvx_precision low
            variables m(p)
            minimize(m'*Sig_eps*m)
            subject to
            l2 <= Sig*m <= u2
        cvx_end
        if( strcmp(cvx_status,'Solved')|strcmp(cvx_status,'Inaccurate/Solved') )
            M(i,:) = m;
            disp('second_solved')
        end
    end
end
M = Q*M;
bet_u = bet_n + M*x'*(y-x*bet_n-int)/n;
V = diag(M*Sig*M');
width = norminv(1-(1-level)/2)*sigma*sqrt(V/n);
CI_u = bet_u + width;
CI_l = bet_u - width;