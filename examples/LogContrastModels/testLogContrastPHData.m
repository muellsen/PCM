%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for soil pH/microbiome data
%
% Both the standard LS and the Huber model as objective function.
% We perform perspective M-estimation for log-contrast regression.
% Model selection is done with theoretical lambda + stability selection
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
       
% Center Y
y_bar = mean(Y);
Y_cent = Y-y_bar;

% Theoretical lambda for model selection
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = sqrt(2/n)*norminv(1-kk/p);

% Linear constraint for log-contrast model
Ceq = ones(1,p); % Standard log-constrast constraint
rhsvec = zeros(size(Ceq,1),1);

% Generate projection matrix that satsifies the constraints
PCeq = Ceq'*pinv(Ceq*Ceq');

%% Perspective estimation with the Huber function

% Algorithmic parameters
clear pcmopts;

% Optimization model
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

pcmopts.abstol = 1e-5;
pcmopts.lamPath = n*lam0;
pcmopts.gamma = 1;

% Model selection with theoretical lambda
t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% Stability selection parameters
ss_perc = 0.5;
n_s = ss_perc*n;
nSS = 100;

% Theoretical lambda for subset selection
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0_ss = sqrt(2/n_s)*norminv(1-kk/p);
pcmopts.lamPath = n_s*lam0_ss;

% Stability selection
[selInds1,selFreq1,selMat1] = stabsel('pcmC2',X,Y_cent,p,5,ss_perc,nSS,pcmopts);


%% Perspective estimation with the LS function

% Optimization model
objFun2 = 'Lq';

pcmopts.objFun = objFun2;
pcmopts.lamPath = n*lam0;

% Model selection with theoretical lambda
t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

pcmopts.lamPath = n_s*lam0_ss;

% Stability selection
[selInds2,selFreq2,selMat2] = stabsel('pcmC2',X,Y_cent,p,5,ss_perc,nSS,pcmopts)


% Plot basic results

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
plot(X*betaPCM1Mat,Y_cent,'.','MarkerSize',10);
title('Huber')
xlabel('Prediction');
ylabel('Data');
grid on

figure;
plot(X*betaPCM2Mat,Y_cent,'.','MarkerSize',10);
title('L2')
xlabel('Prediction');
ylabel('Data');
grid on

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
set(gca,'XTickLabel',jointSet,'fontsize',18)
box on
xtickangle(90)
xlim([0 setLen+0.5])





