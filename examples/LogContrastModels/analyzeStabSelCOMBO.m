%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analysis of log-contrast modeling for prediction of BMI from gut
% microbiome + covariates with standard linear constraints
%
% Both the standard LS and the Huber model are used as objective function.
% The script analyzes the stability selection solution and 
% reproduces Figure 4 in Section 4.1 in [3] using the data from 
% testLogContrastCOMBO.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Two selected sets
selSet1 = find(selFreq1>0.1);
selSet2 = find(selFreq2>0.1);

% Joint sets
jointSet = union(selSet1,selSet2);
setLen = length(jointSet);

% The set of genera is available in PredLabels

% Figure to show the selection probabilities of all selected taxa
figure;stem((1:setLen)-0.25,selFreq2(jointSet),'LineWidth',5);
hold on;stem((1:setLen)+0.25,selFreq1(jointSet),'LineWidth',5);
legend('LS','Huber','Location','NorthWest')
grid on
xlabel('Taxon i')
ylabel('Selection frequency')
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen))
set(gca,'XTickLabel',PredLabels(jointSet),'fontsize',18)
box on
xtickangle(90)
xlim([0.25 setLen+0.5])

% Threshold 
th_stab = 0.7;

selTop1 = find(selFreq1>=th_stab) 
selTop2 = find(selFreq2>=th_stab) 

% Reoptimize over the subset
X_sel = X(:,selTop1);

[n,p] = size(X_sel);

% Center X and Y
y_bar = mean(Y);
Y_cent = Y-y_bar;

% Theoretical lambda (recalculated)
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = sqrt(2/n)*norminv(1-kk/p);


% Linear constraint for log-contrast model
Ceq = ones(1,p); % Standard log-constrast constraint
rhsvec = zeros(size(Ceq,1),1);

% Generate betaTrue that satsifies the constraints
PCeq = Ceq'*pinv(Ceq*Ceq');


% Algorithmic parameters
clear pcmopts;

% Set-up for refitting
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
pcmopts.lamPath = 0;%n*lam0;
pcmopts.gamma = 1;

t1=now;
[betaPCM1Matstab, sigmaPCM1Matstab,funPCM1Mat,outPCM1stab] = pcmC2(X_sel, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% Contraints + no robustness
objFun2 = 'Lq';

pcmopts.objFun = objFun2;

t1=now;
[betaPCM2Matstab, sigmaPCM2Matstab,funPCM2Matstab,outPCM2stab] = pcmC2(X_sel, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

% Inspect taxonomic names at the lowest possible level 
topTaxa = PredLabels(selTop1)

% Insert by hand
topNames = topTaxa

% Figure to show the beta values for the selected robust variables
figure;
stem((1:p),betaPCM2Matstab,'LineWidth',20);
%legend('Huber','Location','NorthWest')
grid on
ylabel('\beta_i','FontSize',40)
set(gca,'FontSize',20)
set(gca,'XTick',(1:p))
set(gca,'XTickLabel',topNames,'fontsize',18)
box on
xtickangle(90)
xlim([0.5 p+0.5])
ylim([-1.1 1.1])

% Identify outliers using the Huber solution
currOut = outPCM1stab;
currBeta = betaPCM1Matstab;
currSigma = sigmaPCM1Matstab;

currInd = 1;

res = (Y_cent-X_sel*currBeta(:,currInd));
eta = currSigma(:,currInd);
[fun,quadInds,linInds] = HuberFun(eta,res,currOut.opts.rho1,currOut.opts.fitLin(1),currOut.opts.qPower1);

% Figure to show the selection probabilities of all selected taxa
figure;
plot(meanBMI+Y(quadInds),meanBMI+X_sel(quadInds,:)*betaPCM1Matstab,'.','MarkerSize',50);
hold on;
plot(meanBMI+Y(linInds),meanBMI+X_sel(linInds,:)*betaPCM1Matstab,'*','MarkerSize',20);
legend('Inliers','Outliers','Location','NorthWest')
grid on
xlabel('BMI measurements','FontSize',40)
ylabel('BMI predictions','FontSize',20)
set(gca,'FontSize',20)
h=line([10,45],[10,45],'Color','black','LineStyle','--','Linewidth',1);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% Summary of the model fit for the Huber model
MAD_all = mean(abs(res));
MAD_in = mean(abs(res(quadInds)));
MAD_out = mean(abs(res(linInds)));

LS_all = sqrt(sum(res.^2)/length(res));
LS_in = sqrt(sum((res(quadInds).^2)/length(quadInds)));
LS_out = sqrt(sum((res(linInds).^2))/length(linInds));

R2 = (corr(Y_cent,X_sel*currBeta(:,currInd))).^2

R2_inliers = (corr(Y_cent(quadInds),X_sel(quadInds,:)*currBeta(:,currInd))).^2
R2_outliers = (corr(Y_cent(linInds),X_sel(linInds,:)*currBeta(:,currInd))).^2


% Identify outliers using the Huber solution
currOut = outPCM2stab;
currBeta = betaPCM2Matstab;
currSigma = sigmaPCM2Matstab;

currInd = 1;

res2 = (Y_cent-X_sel*currBeta(:,currInd));
eta2 = currSigma(:,currInd);
[fun2,quadInds2,linInds2] = HuberFun(eta2,res2,currOut.opts.rho1,currOut.opts.fitLin(1),currOut.opts.qPower1);

% Summary of the model fit for the LS model

MAD_all2 = mean(abs(res2));
MAD_in2 = mean(abs(res2(quadInds2)));
MAD_out2 = mean(abs(res2(linInds2)));

LS_all2 = sqrt(sum(res2.^2)/length(res2));
LS_in2 = sqrt(sum((res2(quadInds2).^2)/length(quadInds2)));
LS_out2 = sqrt(sum((res2(linInds2).^2))/length(linInds2));

R2_2 = (corr(Y_cent,X_sel*currBeta(:,currInd))).^2



