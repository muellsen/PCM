function bet = bycvx(y, x, constr2, lam)
% Use CVX to solve sparse regression (with a linear constraint)

[n p] = size(x);

if (sum(abs(constr2))==0)
	cvx_begin quiet
    cvx_precision low
    variable bet(p)
        minimize( (y-x*bet)'*(y-x*bet)+2*n*lam*norm(bet,1) )
	cvx_end
else
	cvx_begin quiet
    cvx_precision low
    variable bet(p)
        minimize( (y-x*bet)'*(y-x*bet)+2*n*lam*norm(bet,1) )
        subject to
        constr2*bet == 0
	cvx_end
end