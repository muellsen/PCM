%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for COMBO data; prediction of BMI from gut
% microbiome + covariates with general subcompositional constraints
%
% Both the standard LS and the Huber model as objective function.
% We perform perspective M-estimation for log-contrast regression.
% Model selection is done with theoretical lambda + stability selection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random number seed
rng(2342)

% load all csv files
load('data/CaloriData.csv')
load('data/FatData.csv')
load('data/BMI.csv')
load('data/CFiltered.mat')
load('data/GeneraFilteredCounts.csv')

Ceq45 = table2array(CFiltered(1:end,2:5));

% load phylogenetic tree
load('data/GeneraPhylo.mat')	

% BMI data (n=96)
Yunorm = BMI;
meanBMI = mean(BMI);
Y = Yunorm-meanBMI;

% Number of samples
n = length(Y);

X_C = CaloriData-mean(CaloriData);
X_F = FatData-mean(FatData);

% Predictor labels
PredLabels = [CFiltered.VarName1;'Calorie';'Fat'];

% Size of the microbiome data (filtered data)
[pG,nG] = size(GeneraFilteredCounts);

% Count data of 45 genera 
Xunorm = GeneraFilteredCounts;% + rand(pG,nG);

% CLR transform data
X0 = clr(Xunorm,1/2)';

X = [X0,X_C,X_F,ones(n,1)];

[n,p] = size(X);
     
n_orig = n;
p_orig = p;

% Center X and Y
y_bar = mean(Y);
Y_cent = Y-y_bar;

% Original X
X_cent = X;

% Theoretical lambda
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = sqrt(2/n)*norminv(1-kk/p);

% Linear constraint for log-contrast model
Ceq = [Ceq45;zeros(3,4)]'; % Subcompositional constraint log-constrast constraint for the compositions
rhsvec = zeros(size(Ceq,1),1);

% Generate betaTrue that satsifies the constraints
PCeq = Ceq'*pinv(Ceq*Ceq');

% Algorithmic parameters
clear pcmopts;

% Optimization model (Huber)
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
pcmopts.regWeights = [ones(p-3,1);[0;0;0]];

pcmopts.abstol = 1e-6;
pcmopts.lamPath = n*lam0;
pcmopts.gamma = 1;

t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

ss_perc = 0.5;
n_s = ss_perc*n;
nSS = 100;

% Theoretical lambda for subset selection
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0_ss = sqrt(2/n_s)*norminv(1-kk/p);
pcmopts.lamPath = n_s*lam0_ss;

[selInds1,selFreq1,selMat1] = stabsel('pcmC2',X_cent,Y_cent,p,5,ss_perc,nSS,pcmopts);

% Optimization model (L2)
objFun2 = 'Lq';

pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;

t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

pcmopts.lamPath = n_s*lam0_ss;

[selInds2,selFreq2,selMat2] = stabsel('pcmC2',X_cent,Y_cent,p,5,ss_perc,nSS,pcmopts);

close all

figure;
stem(betaPCM1Mat,'LineWidth',1);
legend('Huber')
xlabel('indices');
ylabel('\beta_i');
grid on

figure;
stem(betaPCM2Mat,'LineWidth',1);
xlabel('indices');
ylabel('\beta_i');
legend('L2')
grid on

figure;
plot(Y_cent+meanBMI,X_cent*betaPCM1Mat+meanBMI,'.','MarkerSize',40);
title('Huber')
xlabel('Data');
ylabel('Prediction');
grid on
xlim([10 max(Y_cent+meanBMI)])
ylim([10 max(Y_cent+meanBMI)])

R2 = (corr(Y_cent,X_cent*betaPCM1Mat)).^2

figure;
plot(Y_cent+meanBMI,X_cent*betaPCM2Mat+meanBMI,'.','MarkerSize',40);
title('L2')
xlabel('Data');
ylabel('Prediction');
grid on
xlim([10 max(Y_cent+meanBMI)])
ylim([10 max(Y_cent+meanBMI)])

R2_2 = (corr(Y_cent,X_cent*betaPCM2Mat)).^2

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
pcmopts.regWeights = [ones(p-3,1);[0;0;0]];
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


% Compare this to the standard constraint without subcompositions

% Linear constraint for log-contrast model
Ceq2 = [ones(1,pG),zeros(1,3)]; % Standard log-constrast constraint for the compositions
rhsvec = zeros(size(Ceq2,1),1);

% Generate betaTrue that satsifies the constraints
PCeq = Ceq2'*pinv(Ceq2*Ceq2');

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
pcmopts.Ceq = Ceq2;
pcmopts.rhsvec = rhsvec;
pcmopts.regWeights = [ones(p-3,1);[0;0;0]];

pcmopts.abstol = 1e-6;
pcmopts.lamPath = n*lam0;
pcmopts.gamma = 1;

t1=now;
[betaPCM1Matones, sigmaPCM1Matones,funPCM1Matones,outPCM1ones] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

ss_perc = 0.5;
n_s = ss_perc*n;
nSS = 100;

% Theoretical lambda for subset selection
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0_ss = sqrt(2/n_s)*norminv(1-kk/p);
pcmopts.lamPath = n_s*lam0_ss;

[selInds1ones,selFreq1ones,selMat1ones] = stabsel('pcmC2',X_cent,Y_cent,p,5,ss_perc,nSS,pcmopts);

% Contraints + no robustness
objFun2 = 'Lq';

pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;

t1=now;
[betaPCM2Matones, sigmaPCM2Matones,funPCM2Matones,outPCM2ones] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

pcmopts.lamPath = n_s*lam0_ss;

[selInds2ones,selFreq2ones,selMat2ones] = stabsel('pcmC2',X_cent,Y_cent,p,5,ss_perc,nSS,pcmopts);



