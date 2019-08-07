% Test functionality of prox solver for different loss functions
close all

% Dimension = number of columns of X
p = 30;
disp(['Dimension of predictors: ',num2str(p)])

% Number of non-zero entries of the regression vector
nnzs=6;
disp(['Number of nnzs: ',num2str(nnzs)])

% Sample sizes
alpha = 3;
n = round(alpha*nnzs*log(p));
disp(['Number of data points: ',num2str(n)])

% Generate leading non-zero entries of size one
const = 1;
firstEntries = const*[(-1).^(1:nnzs)]';

% True beta vector
betaTrue = [firstEntries;zeros(p-nnzs,1)]; % only nnz nonzero coefficients

% Correlation kappa
kappa = 0.7;

% Generate covariance matrix
covMat = kappa*ones(p,p);
covMat(1:p+1:p^2) = 1;
cholCov = chol(covMat);

%  Noise vector
n_1 = round(n/4);
n_2 = round(n/4);
n_3 = n-n_1-n_2;

s1 = 5;
s2 = 0.5;
s3 = 0.01;

sig1 = s1*ones(n_1,1);
sig2 = s2*ones(n_2,1);
sig3 = s3*ones(n_3,1);

sig = ones(n,1);

% Generate data
X = (cholCov'*randn(p,n))';

% Normalize X to length sqrt(n)
normX = repmat(sqrt(sum(X.^2)),n,1);
X = sqrt(n)*X./normX;

% Gaussian noise
noiseVec = randn(n,1);

% Response with sigma * standardized noise vector
Y = X*betaTrue + [sig1;sig2;sig3].*noiseVec;
%Y = X*betaTrue + sig.*noiseVec;

% Outliers
nOutl = n/4;
outlInds = unique(randi(n,round(nOutl),1));
%Y(outlInds) = max(Y)*randn(length(outlInds),1);

t1=now;
[betaSCMat, sigmaSCMat,funSCMat,outSC] = sc_lasso(X, Y);
t2=now;
timeSCLasso = (t2-t1)*(60*60*24)

% Square case
qPower1 = 2;

% Algorithmic parameters
clear pcmopts;
pcmopts.abstol = 1e-4;
pcmopts.lambda = n*outSC.lamPath(1:10:end);

pcmopts.objFun = 'Huber';
pcmopts.qPower1 = 2;
pcmopts.fitLin = 1/2;
pcmopts.nGroupSiz = [n_1,n_2,n_3];

pcmopts.penFun = 'BerHu';
pcmopts.qPower2 = 2;
pcmopts.penLin = 1/2;

pcmopts.regFun = '';
pcmopts.gamma = 1;

t1=now;
[betaPCMMat,sigmaPCMMat,tauPCMMat, outPCM] = pcmC2(X, Y, pcmopts);
t2=now;

timePCM = (t2-t1)*(60*60*24)

clear huberopts;
huberopts.rho = 1.345;  % Huber parameter
huberopts.lambda = 2*n*outSC.lamPath(1:10:end);
%huberopts.lenPath = lenPath; % Path length

t1=now;
[betaHuberMat, sigmaHuberMat,funHuberMat,outHuber] = huber_cvx(X, Y, huberopts);
t2=now;
timeHuber = (t2-t1)*(60*60*24)


% Algorithmic parameters
clear pcmopts;
pcmopts.abstol = 1e-4;
pcmopts.lambda = n*outSC.lamPath(1:10:end);

pcmopts.objFun = 'Lq';
pcmopts.qPower1 = 2;
pcmopts.fitLin = 1/2;
pcmopts.nGroupSiz = [n_1,n_2,n_3];

pcmopts.penFun = '';
pcmopts.qPower2 = 2;
pcmopts.penLin = 1/2;

pcmopts.regFun = 'L1';
pcmopts.gamma = 1;

t1=now;
[betaPCMMat2,sigmaPCMMat2,tauPCMMat2, outPCM2] = pcmC2(X, Y, pcmopts);
t2=now;

timePCM2 = (t2-t1)*(60*60*24)

% Algorithmic parameters
clear pcmopts;
pcmopts.abstol = 1e-4;
pcmopts.lambda = n*outSC.lamPath(1:10:end);

pcmopts.objFun = 'Huber';
pcmopts.qPower1 = 2;
pcmopts.fitLin = 1/2;
%pcmopts.nGroupSiz = [n_1,n_2,n_3];

pcmopts.penFun = '';
pcmopts.qPower2 = 2;
pcmopts.penLin = 1/2;

pcmopts.regFun = 'L1';
pcmopts.gamma = 1;

t1=now;
[betaPCMMat3,sigmaPCMMat3,tauPCMMat3, outPCM3] = pcmC2(X, Y, pcmopts);
t2=now;

timePCM3 = (t2-t1)*(60*60*24)

figure;
plot(outSC.lamPath,betaSCMat')
title('Scaled Lasso')

figure;
plot(outPCM.lamPath,betaPCMMat')
title('Het-Huber + BerHu')

figure;
plot(outPCM2.lamPath,betaPCMMat2')
title('Het-Lasso + L1')

figure;
plot(outHuber.lamPath,betaHuberMat')
title('Huber + L1(cvx)')

figure;
plot(outPCM3.lamPath,betaPCMMat3')
title('Huber + L1')

figure;
plot(outPCM2.lamPath,sigmaPCMMat2','LineWidth',4)
hold on
plot(outPCM2.lamPath,s1*ones(1,length(outPCM2.lamPath)),'k--','LineWidth',4)
plot(outPCM2.lamPath,s2*ones(1,length(outPCM2.lamPath)),'k--','LineWidth',4)
plot(outPCM2.lamPath,s3*ones(1,length(outPCM2.lamPath)),'k--','LineWidth',4)
ylim([0 s1+1])
grid on
title('Heteroscedastic Scaled Lasso + L1')

figure;
plot(outPCM.lamPath,2*sigmaPCMMat','LineWidth',4)
hold on
plot(outPCM.lamPath,s1*ones(1,length(outPCM.lamPath)),'k--','LineWidth',4)
plot(outPCM.lamPath,s2*ones(1,length(outPCM.lamPath)),'k--','LineWidth',4)
plot(outPCM.lamPath,s3*ones(1,length(outPCM.lamPath)),'k--','LineWidth',4)
ylim([0 s1+1])
grid on
title('Heteroscedastic Huber + BerHu')








