%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for COMBO data; prediction of BMI from gut
% microbiome + covariates with standard linear constraints
%
% Both the standard LS and the Huber model as objective function.
% We perform perspective M-estimation for log-contrast regression.
% This script computes  
% 1) Stability selection with theoretical lambda 
% 2) Regularization path on all data + theoretical lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random number seed
rng(2342)

% load all csv files
load('data/COMBO/CaloriData.csv')
load('data/COMBO/FatData.csv')
load('data/COMBO/BMI.csv')
load('data/COMBO/GeneraFilteredCounts.csv')
load('data/COMBO/GeneraCounts.csv')
load('data/COMBO/CFiltered.mat')

Ceq45 = table2array(CFiltered(1:end,2:5));

% load phylogenetic tree
load('data/COMBO/GeneraPhylo.mat')

% BMI data (n=96)
Yunorm = BMI;
meanBMI = mean(BMI);
Y = Yunorm-meanBMI;

% Number of samples
n = length(Y);

% Covariate data
X_C = CaloriData-mean(CaloriData);
X_F = FatData-mean(FatData);

% Predictor labels
PredLabels = [GeneraPhylo(:,7);'Calorie';'Fat'];

% Size of the microbiome data
[pG,nG] = size(GeneraCounts);

% Count data of 87 genera 
Xunorm = GeneraCounts;

% CLR transform data with pseudo count of 0.5
X0 = clr(Xunorm,1/2)';

% Joint microbiome and covariate data and offset
X = [X0,X_C,X_F,ones(n,1)];

% New dimensionality
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

% Linear constraint for log-contrast model (on microbiome data)
Ceq = [ones(1,p-3),0,0,0]; % Standard log-constrast constraint for the compositions
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

% Optimization model (LS)

objFun2 = 'Lq';

pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;

t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

pcmopts.lamPath = n_s*lam0_ss;

[selInds2,selFreq2,selMat2] = stabsel('pcmC2',X_cent,Y_cent,p,5,ss_perc,nSS,pcmopts);

% Plot summary of results so far
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
xlabel('Prediction');
ylabel('Data');
grid on
xlim([10 max(Y_cent+meanBMI)])
ylim([10 max(Y_cent+meanBMI)])

R2 = (corr(Y_cent,X_cent*betaPCM1Mat)).^2

figure;
plot(Y_cent+meanBMI,X_cent*betaPCM2Mat+meanBMI,'.','MarkerSize',40);
title('L2')
xlabel('Prediction');
ylabel('Data');
grid on
xlim([10 max(Y_cent+meanBMI)])
ylim([10 max(Y_cent+meanBMI)])

R2_2 = (corr(Y_cent,X_cent*betaPCM2Mat)).^2

% Selection frequency
selSet1 =find(selFreq1~=0);
selSet2 =find(selFreq2~=0);

jointSet = union(selSet1,selSet2);
setLen = length(jointSet);

figure;stem((1:setLen)-0.5,selFreq2(jointSet),'LineWidth',5);hold on;stem(1:setLen,selFreq1(jointSet),'LineWidth',5);
legend('LS','Huber','Location','NorthWest')
grid on
xlabel('Taxon i')
ylabel('Selection frequency')
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen)-0.25)
set(gca,'XTickLabel',PredLabels(jointSet),'fontsize',18)
box on
xtickangle(90)
xlim([0 setLen+0.5])


% Compute the entire solution path on all data for reference
clear pcmopts

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
pcmopts.regWeights = [ones(p-3,1);[0;0;0]]; % do not regularize covariate and offset
pcmopts.lenPath = 40;

pcmopts.abstol = 1e-6;
pcmopts.gamma = 1;

t1=now;
[betaPCM1Matpath, sigmaPCM1Matpath,funPCM1Matpath,outPCM1path] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% % Optimization model (L2)
objFun2 = 'Lq';
pcmopts.objFun = objFun2;

t1=now;
[betaPCM2Matpath, sigmaPCM2Matpath,funPCM2Matpath,outPCM2path] = pcmC2(X_cent, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)
