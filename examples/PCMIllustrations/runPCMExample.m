% Run PCM solver for different loss function on the same data set with
% extension to outlier removal

rng('default')

rng(23)

close all

% Dimension = number of columns of X
p = 64;
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
kappa = 0.3;

% Generate covariance matrix
covMat = kappa*ones(p,p);
covMat(1:p+1:p^2) = 1;
cholCov = chol(covMat);

%  Noise vector
n_1 = round(n/3);
n_2 = round(n/3);
n_3 = n-n_1-n_2;

s1 = 5;
s2 = 0.5;
s3 = 0.1;

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
Y1 = X*betaTrue + [sig1;sig2;sig3].*noiseVec;
%Y = X*betaTrue + sig.*noiseVec;

% Outliers
nOutl = round(0.1*n);
temp = randperm(n);
outlInds = temp(1:nOutl);
oVec = zeros(n,1);
oVec(outlInds,1) = mean(Y1)+5*randn(nOutl,1);

Y = Y1 + oVec;

figure;plot(1:n,Y,'.','MarkerSize',20);hold on;plot(outlInds,Y(outlInds),'.','MarkerSize',20)
figure;plot(Y,X*betaTrue,'.','MarkerSize',20);hold on;plot(Y(outlInds),X(outlInds,:)*betaTrue,'.','MarkerSize',20)

% Baseline method (Coordinate descent)
scopts.lenPath = 50;
t1=now;
[betaSCMat, sigmaSCMat,funSCMat,outSC] = sc_lasso(X, Y,scopts);
t2=now;
timeSCLasso = (t2-t1)*(60*60*24)

% Set=up experimental design
power1Vec = [2,3/2];
nPow1 = length(power1Vec);

objFunCell = {'Lq','Huber'}
nObj = length(objFunCell);

penFunCell = {'','Lq','BerHu'}
nPen = length(penFunCell);

power2Vec = [2];
nPow2 = length(power2Vec);

regCell = {'','L1','L1s'}
nReg = length(regCell);

modes = {'homosced','hetsced'}
nModes = length(modes);

nExperiments = nObj*nPow1*nPen*nPow2*nReg*nModes;

betaCell = cell(nExperiments,1);
sigmaCell = cell(nExperiments,1);
tauCell = cell(nExperiments,1);
outCell = cell(nExperiments,1);
runTimeVec = zeros(nExperiments,1);

cnt = 1;

for i1 = 1:nPow1
    for i2=1:nObj
        for i3 = 1:nPen
            for i4 = 1:nPow2
                for i5 = 1:nReg
                    for i6 = 1:nModes
                        
                        % Algorithmic parameters reset
                        clear pcmopts;
                        
                        % New set of parameters
                        pcmopts.abstol = 5e-4;
                        
                        if strcmp(penFunCell{i3},'BerHu')
                            pcmopts.lamPath = n/2*outSC.lamPath;
                        else
                            pcmopts.lamPath = n*outSC.lamPath;
                        end
                        
                        pcmopts.qPower1 = power1Vec(i1);
                        
                        pcmopts.objFun = objFunCell{i2};
                        pcmopts.rho1 = 1.345;
                        pcmopts.fitLin = 1/2;
                        if strcmp(modes{i6},'hetsced')
                            pcmopts.nGroupSiz = [n_1,n_2,n_3];
                        end
                        
                        pcmopts.penFun = penFunCell{i3};
                        pcmopts.rho2 = 1.345;
                        pcmopts.qPower2 = power2Vec(i4);
                        pcmopts.penLin = 1/2;
                        %pcmopts.penLB = 0.5;
                        
                        pcmopts.regFun = regCell{i5};
                        pcmopts.gamma = 1
                        
                        t1=now;
                        [betaPCMMat,sigmaPCMMat,tauPCMMat, outPCM] = pcmC2(X, Y, pcmopts);
                        t2=now;
                        timePCM = (t2-t1)*(60*60*24)
                        
                        betaCell{cnt} = betaPCMMat;
                        sigmaCell{cnt} = sigmaPCMMat;
                        tauCell{cnt} = tauPCMMat;
                        outCell{cnt} = outPCM;
                        runTimeVec(cnt) = timePCM;
                        
                        disp([num2str(cnt/nExperiments*100),' % percent done'])
                        
                        if strcmp(penFunCell{i3},'') && strcmp(regCell{i5},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1))];
                            
                        elseif strcmp(penFunCell{i3},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Reg: ',regCell{i5},', Mode:', modes{i6}];
                            
                        elseif strcmp(regCell{i5},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Pen: ' penFunCell{i3}, ' with q=',num2str(power2Vec(i4)), ', Mode:', modes{i6}];
                        else
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Pen: ' penFunCell{i3}, ' with q=',num2str(power2Vec(i4)), ', Reg: ',regCell{i5},', Mode:', modes{i6}];
                        end
                        
                        currBeta = betaCell{cnt};
                        currSigma = sigmaCell{cnt};
                        currTau = tauCell{cnt};
                        currOut = outCell{cnt};
                        currRuntime = runTimeVec(cnt);
                        
                        figure(cnt);
                        semilogx(currOut.lamPath,currBeta,'LineWidth',1)
                        hold on
                        semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        xlabel('Regularization parameter \lambda')
                        ylabel('Solution path')
                        grid on
                        title(titleString)
                        drawnow;
                        
                        cnt = cnt+1;
                        
                        save('IllustEx_p64_2')
                    end
                end
            end
        end
    end
end










