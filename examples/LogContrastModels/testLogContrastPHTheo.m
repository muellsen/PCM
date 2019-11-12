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
% This influences the results slightly 
% due to the randomnes in stability selection

rng(2342)

% Load pH data Y and microbiome data X (already log-transformed)
load('data/pHData.mat')
[n,p] = size(X);

% Match OTUs from Morton et al. and own analysis
matchOTU_IDS;
       
% Center X and Y
y_bar = mean(Y);
Y_cent = Y-y_bar;

loc_x = mean(X);
sca_X = std(X);
b2 = 1./sca_X;
diagB = diag(b2);
X_c = (X - ones(n,1)*loc_x)*diagB;
X_cent = X_c;

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
pcmopts.gamma = 1;

t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% Contraints + no robustness
objFun2 = 'Lq';

pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;

t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X_cent, Y_cent, pcmopts);
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
pcmopts.gamma = 1;

t1=now;
[betaPCM1Matpath, sigmaPCM1Matpath,funPCM1Matpath,outPCM1path] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% L2 
objFun2 = 'Lq';

pcmopts.objFun = objFun2;

t1=now;
[betaPCM2Matpath, sigmaPCM2Matpath,funPCM2Matpath,outPCM2path] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

