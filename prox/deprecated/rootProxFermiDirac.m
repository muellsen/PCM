% Testing how to find the root for Fermi-Dirac

gamma = 2;
eta = 1;

y=randn;

pFD = @(p) gamma*p + (eta+log(1+exp(p)))/(1+exp(-p))-y;

tic
p = fzero(pFD,0,[-1,1])
toc


