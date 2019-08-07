function [ BICVec ] = BIC_Huber(X,Y,sigmaVec,betaMat,rho)
%BIC_Huber BIC for Huber-type penalized regression
%   Detailed explanation goes here
%
% Dimension of the problem
[n,p] = size(X)

% Length of regularization path
nLams = length(sigmaVec);

% Thresholding for numerical zeros
epsNum = 1e-6;

% # of non-zero elements
kVec = sum(abs(betaMat)>epsNum)

% DoG term in BIC 
degTerm = kVec.*log(n)/(2*n);

% 
BICVec = zeros(nLams,1);

for i=1:nLams
    
    tempHuber=0;
    for j=1:n
        tempHuber = tempHuber+objHuber1D(X(j,:),Y(j),betaMat(:,i),sigmaVec(i),rho);
    end
    
    BICVec(i) = tempHuber;
end

% Full Huber BIC
BICVec = BICVec + degTerm;

function huberFun = objHuber1D(X,Y,betaVec,rho)

res = (Y-X*betaVec)









