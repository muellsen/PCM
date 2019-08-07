% Test functionality of prox solver

% Dimension = number of columns of X
p = 50;
disp(['Dimension of predictors: ',num2str(p)])

% Sample sizes
n = 20;
disp(['Number of data points: ',num2str(n)])

% Number of non-zero entries of the regression vector
nnzs=5;
disp(['Number of nnzs : ',num2str(nnzs)])

% Generate leading non-zero entries of size one
firstEntries = 1*[(-1).^(1:nnzs)]';

% True beta vector
betaTrue = [firstEntries;zeros(p-nnzs,1)]; % only nnz nonzero coefficients

% Correlation kappa
kappa = 0;

% Generate covariance matrix
covMat = kappa*ones(p,p);
covMat(1:p+1:p^2) = 1;
cholCov = chol(covMat);

%  noise vector
sig = 1;

% Generate data
X = (cholCov'*randn(p,n))';

% Normalize X to length sqrt(n)
normX = repmat(sqrt(sum(X.^2)),n,1);
X = sqrt(n)*X./normX;

% Gaussian noise
noiseVec = randn(n,1);

% Response with sigma * standardized noise vector
Y = X*betaTrue + sig*noiseVec;

% Regularization path for c
cpath = 0.5;%[0.7:-0.05:0.05];


cntError = 0;

nTests = 1e3;
proxDiffVec = zeros(nTests,1)

for i=1:nTests
    
    % random beta
    betaTest1 = randn(p,1);
    betaTest2 = randn(p,1);
    
    Yvec = Y;
    iRand = randi(p);  
    v = X(:,iRand);
    eta1 = v'*X*betaTest1;
%     
%     if eta1<0
%         v = -v;
%         eta1 = v'*X*betaTest1;
%     end
    y1 = X*betaTest1;

    %eta1 = rand;
    %y1 = randn(n,1);
    
    [etaProx1,yProx1] = proxgTREX(eta1,y1,v,Yvec);
    %[etaProx1,yProx1] = proxg0TREX(eta1,y1);
    
    
    eta2 = v'*X*betaTest2;
   
%     if eta2<0
%         v = -v;
%         eta2 = v'*X*betaTest2;
%     end
%    
    y2 = X*betaTest2;
    
    
    %eta2 = rand;
    %y2 = randn(n,1);
    
    [etaProx2,yProx2] = proxgTREX(eta2,y2,v,Yvec);
    %[etaProx2,yProx2] = proxg0TREX(eta2,y2);

    
    normY = sum(([eta1;y1]-[eta2;y2]).^2);
    
    normDiffY = sum((([eta1;y1]-[etaProx1;yProx1]) - ([eta2;y2]-[etaProx2;yProx2])).^2);
    
    normProxY = sum(([etaProx1;yProx1]-[etaProx2;yProx2]).^2);
    
    
    proxDiff = normY-(normDiffY+normProxY)
    
    proxDiffVec(i,1) = proxDiff;
    
    figure(1);
    plot(1:nTests,proxDiffVec,'LineWidth',5)
    grid on
    xlabel('Iteration')
    ylabel('Diff')
    drawnow
    
    if proxDiff < 0
        
        normY
        normDiffY
        normProxY
        cntError = cntError + 1;
        warning('Formula not satisifed')
        pause
    end
end












