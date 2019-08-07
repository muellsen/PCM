% Run PCM solver for low dimensional example loss function on the same data set with
% extension to outlier removal

rng('default')

rng(23)

close all

% Dimension = number of columns of X
p = 3;
disp(['Dimension of predictors: ',num2str(p)])

% Number of non-zero entries of the regression vector
nnzs=2;
disp(['Number of nnzs: ',num2str(nnzs)])

% Sample sizes
alpha = 8;
n = round(alpha*nnzs*log(p));
disp(['Number of data points: ',num2str(n)])

% Generate leading non-zero entries of size one
const = 0.25;
firstEntries = const*[(-1).^(1:nnzs)]';

% True beta vector
betaTrue = [firstEntries;zeros(p-nnzs,1)]; % only nnz nonzero coefficients

% Generate covariance matrix
sig = 3*ones(n,1);

% Correlation kappa
kappa = 0;

% Generate covariance matrix
covMat = kappa*ones(p,p);
covMat(1:p+1:p^2) = 1;
cholCov = chol(covMat);

% Generate data
X = (cholCov'*randn(p,n))';

% Normalize X to length sqrt(n)
normX = repmat(sqrt(sum(X.^2)),n,1);
X = sqrt(n)*X./normX;

n_1 = round(n/2);
n_2 = n-n_1;

% Gaussian noise
noiseVec = [randn(n_1,1);zeros(n_2,1)];

oVec = zeros(n,1);
oVec(end) = 1;

% Response with sigma * standardized noise vector
Y1 = X*betaTrue+sig.*noiseVec; 
Y2 = Y1 + oVec;


% Baseline method (Coordinate descent)
scopts.lenPath = 200;
t1=now;
[betaSCMat, sigmaSCMat,funSCMat,outSC] = sc_lasso(X, Y2,scopts);
t2=now;
timeSCLasso = (t2-t1)*(60*60*24)

figure;plot(outSC.lamPath,sigmaSCMat)
figure;plot(outSC.lamPath,betaSCMat')
drawnow

% Set=up experimental design
power1Vec = [2];
nPow1 = length(power1Vec);

objFunCell = {'Lq','Huber'}
nObj = length(objFunCell);

lbVec = [0,0.05];
nLB = length(lbVec);

penFunCell = {''};

dataType = [Y1,Y2];
nData = size(dataType,2);

regCell = {'L1'}
nReg = length(regCell);

modes = {'homosced','hetsced'}
nModes = length(modes);

nExperiments = nObj*nPow1*nLB*nData*nReg*nModes;

betaCell = cell(nExperiments,1);
sigmaCell = cell(nExperiments,1);
tauCell = cell(nExperiments,1);
outCell = cell(nExperiments,1);
runTimeVec = zeros(nExperiments,1);

cnt = 1;

for i1 = 1:nPow1
    for i2=1:nObj
        for i3 = 1:nLB
            for i4 = 1:nData
                for i5 = 1:nReg
                    for i6 = 1:nModes
                        
                        % Algorithmic parameters reset
                        clear pcmopts;
                        
                        % New set of parameters
                        pcmopts.abstol = 5e-7;
                        pcmopts.lamPath = sqrt(3)*n*outSC.lamPath;
                        
                        pcmopts.fitLB = lbVec(i3);
                        
                        pcmopts.qPower1 = power1Vec(i1);
                        
                        pcmopts.objFun = objFunCell{i2};
                        pcmopts.rho1 = 1.345;
                        pcmopts.fitLin = 1/2;
                        if strcmp(modes{i6},'hetsced')
                            pcmopts.nGroupSiz = [n_1,n_2];
                        end
                        
                        pcmopts.penFun = '';
                        
                        pcmopts.regFun = regCell{i5};
                        pcmopts.gamma = 1;
                        
                        % Data with or without outliers
                        Y = dataType(:,i4);
                        
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
                        
                        % With or without outlier
                        if i4==2
                            outlString = 'Yes';
                        else
                            outlString = 'No';
                        end
                            
                        
                        titleString = ['Fit: ', objFunCell{i2},', Reg: ',regCell{i5},', Outlier: ',outlString,', Mode: ', modes{i6},', LB: ', num2str(lbVec(i3))];
                        
                        currBeta = betaCell{cnt};
                        currSigma = sigmaCell{cnt};
                        currTau = tauCell{cnt};
                        currOut = outCell{cnt};
                        currRuntime = runTimeVec(cnt);
                        
                        figure(cnt);
                        plot(currOut.lamPath,currBeta,'LineWidth',5)
                        hold on
                        %semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        xlabel('Regularization parameter \alpha_1','FontSize',30)
                        ylabel('\beta_i','FontSize',30)
                        grid on
                        title(titleString)
                        drawnow;
                        
                        figure(cnt+23);
                        plot(currOut.lamPath,currSigma,'LineWidth',5)
                        hold on
                        %semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        xlabel('Regularization parameter \alpha_1','FontSize',30)
                        ylabel('\sigma_i','FontSize',30)
                        grid on
                        title(titleString)
                        drawnow;
                        
                        cnt = cnt+1;
                        
                        %save('NoiseFreeOutlier')
                    end
                end
            end
        end
    end
end





