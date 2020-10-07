%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for soil pH/microbiome data
%
% Both the standard LS and the Huber model as objective function.
% We perform perspective M-estimation for log-contrast regression.
% Model selection is done with theoretical lambda alone (Appendix D in [3])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random number seed (for reproducibility)
rng(2342)

% Load pH data Y and microbiome count data X_count
load('data/pHSoilData/pHSoilData.mat')
% load('data/pHSoilData/deprecated/pHDataOld.mat')

% X_count comprises the predictor matrix as [nxp] matrix
[n,p] = size(X_count);

% Count data of 116 genera transformed into pxn matrix
Xunorm = X_count';

% CLR transform data with pseudo count of 0.5
X = clr(Xunorm,1/2)';
       
% Center Y
y_bar = mean(Y);
Y_cent = Y-y_bar;

% Theoretical lambda
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = sqrt(2/n)*norminv(1-kk/p);

% Linear constraint for log-contrast model
Ceq = [ones(1,p)]; % Standard log-constrast constraint for the compositions
rhsvec = zeros(size(Ceq,1),1);

% Generate betaTrue that satsifies the constraints
PCeq = Ceq'*pinv(Ceq*Ceq');

% Algorithmic parameters
clear pcmopts;

% Set-up experimental design
power1 = 2;
objFun = 'Huber';
penFun = 'L1';
regFun = 'Ceq';

pcmopts.qPower1 = power1;
pcmopts.objFun = objFun;
pcmopts.fitLin = 1/2;
pcmopts.penFun = penFun;
pcmopts.regFun = regFun;
pcmopts.Ceq = Ceq;
pcmopts.rhsvec = rhsvec;

pcmopts.abstol = 1e-6;
pcmopts.lamPath = n*lam0;
pcmopts.gamma = 0.1;

t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)


% Contraints + no robustness
objFun2 = 'Lq';
pcmopts.fitLin = 1/2;
pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;
pcmopts.gamma = 0.6;
pcmopts.verbose = 1;

t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

% Compute the entire solution path on all data for reference
clear pcmopts

% Huber
power1 = 2;
objFun = 'Huber';
penFun = 'L1';
regFun = 'Ceq';

pcmopts.qPower1 = power1;
pcmopts.objFun = objFun;
pcmopts.fitLin = 1/2;
pcmopts.penFun = penFun;
pcmopts.regFun = regFun;
pcmopts.Ceq = Ceq;
pcmopts.rhsvec = rhsvec;
pcmopts.lenPath = 40;

pcmopts.abstol = 1e-6;
pcmopts.gamma = 0.5;
pcmopts.verbose =1

t1=now;
[betaPCM1Matpath, sigmaPCM1Matpath,funPCM1Matpath,outPCM1path] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM1path = (t2-t1)*(60*60*24)

% L2 
objFun2 = 'Lq';
pcmopts.objFun = objFun2;
pcmopts.gamma = 2;
pcmopts.verbose =1

t1=now;
[betaPCM2Matpath, sigmaPCM2Matpath,funPCM2Matpath,outPCM2path] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM2path = (t2-t1)*(60*60*24)

