function x = gg6(n,N,mu,beta,rho)
% Generate n Generalized Gaussian m-dimensional vectors with means mu(:),
% inverse scales beta(:), and parameters rho(:).
if length(mu) == 1
    mu = mu*ones(1,n);
    beta = beta*ones(1,n);
    rho = rho*ones(1,n);
end
x = zeros(n,N);
for i = 1:n
   x(i,:) = mu(i) + (1/sqrt(beta(i))) * ( gamrnd(1/rho(i),1,1,N)).^(1/rho(i) ) .* ((rand(1,N)<0.5)*2-1) ;
end   