function bet = lasso_constr(y, x, constr2, lam, tol, maxiter,varargin)

    
% use coordinate descent to solve regularized regression
[n p] = size(x);
k = size(constr2, 1);

gramC = constr2'*constr2;
gramX = x'*x;
diagC = diag(gramC);
diagX = diag(gramX);
dediagC = gramC - diag(diagC);
dediagX = gramX - diag(diagX);
covXY = x'*y;

mu = 1;
bet = ones(p,1)/p;

if isempty(varargin)
    bet0 = zeros(p,1);
else
    bet0 = cell2mat(varargin{1});
end
iter = 0;

if (sum(abs(constr2))==0) % no constraint
    term0 = (covXY - dediagX*bet)/n;
    term2 = diagX/n;
    while (sum(abs(bet-bet0))>tol & iter<maxiter)
        bet0 = bet;
        for j=1:p
            term1 = sign(term0(j)) * max(0, abs(term0(j))-lam);
            bet(j) = term1/term2(j);
            dif = bet(j)-bet0(j);
            term0 = term0 - dediagX(:,j)*dif/n;
        end        
        iter = iter + 1;
    end
    if(iter<maxiter) 
        disp('solved') 
    end
else % with constraint
    ksi = zeros(k,1);
    ksi0 = ones(k,1);
    term0 = (covXY - dediagX*bet)/n - mu*(constr2'*ksi + dediagC*bet);
    term2 = diagX/n + mu*diagC;
    while (sum(abs(ksi-ksi0))>tol && iter<maxiter)
        ksi0 = ksi;
        iter2 = 0;
        bet0 = bet0 + 1;
        while (sum(abs(bet-bet0))>tol && iter2<1000)
            bet0 = bet;
            for j=1:p
                term1 = sign(term0(j)) * max(0, abs(term0(j))-lam);
                bet(j) = term1/term2(j);
                dif = bet(j)-bet0(j);
                term0 = term0 - dediagX(:,j)*dif/n - dediagC(:,j)*dif*mu;
            end
            iter2 = iter2 + 1;
        end
        dif2 = constr2*bet;
        ksi = ksi + dif2; 
        term0 = term0 - mu*constr2'*dif2;
        iter = iter + 1;
    end
    if(iter<maxiter) 
        disp('solved') 
    end
end
