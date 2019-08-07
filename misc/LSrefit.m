function [betaMatRF,L1Vec,L2Vec,R2Vec,nnzVec] = LSrefit(betaMat,X,Y)

[n,p] = size(X);
nLams = size(betaMat,2);

betaMatRF = zeros(p,nLams);
L1Vec = zeros(1,nLams);
L2Vec = zeros(1,nLams);
R2Vec = zeros(1,nLams);
nnzVec = zeros(1,nLams);

for i = 1:nLams
    
    inds = find(betaMat(:,i)~=0);
    S_samp = X(:,inds)'*X(:,inds);
    if cond(S_samp)<1e12
        
        temp = inv(S_samp)*X(:,inds)'*Y;
        betaMatRF(inds,i) = temp;
        L1Vec(i) = sum(abs(Y-X*betaMatRF(:,i)))/n;
        L2Vec(i) = sum((Y-X*betaMatRF(:,i)).^2)/n;
        R2Vec(i) = (corr(Y,X*betaMatRF(:,i))).^2;
    end
    nnzVec(i) = length(inds);
end


