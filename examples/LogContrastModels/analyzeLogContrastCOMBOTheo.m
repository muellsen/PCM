%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analysis of log-contrast modeling for prediction of BMI from gut
% microbiome + covariates with standard linear constraints
%
% Both the standard LS and the Huber model are used as objective function.
% The script analyzes the theoretical solution and 
% reproduces Figure 3 in Section 4.1 in [3] using the data from 
% testLogContrastCOMBO.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thresholding the numerical solutions
epsThresh = 1e-3;
betaPCM1Matpath(abs(betaPCM1Matpath(:))<epsThresh)=0;
betaPCM2Matpath(abs(betaPCM2Matpath(:))<epsThresh)=0;

% Plot the subset of the relevant regularization path
figure;
plot(outPCM1path.lamPath(2:25)/n,betaPCM1Matpath(1:p-3,2:25),'LineWidth',3)
grid on
box on
xlabel('\lambda path')
ylabel('\beta_i')
ylim([-1.5 1.5])
set(gca,'FontSize',20)
hold on
plot(lam0,betaPCM1Mat(1:p-3),'k*','MarkerSize',5)

figure;
plot(outPCM2path.lamPath(2:25)/n,betaPCM2Matpath(1:p-3,2:25),'LineWidth',3)
grid on
box on
xlabel('\lambda path')
ylabel('\beta_i')
ylim([-1.5 1.5])
set(gca,'FontSize',20)
hold on
plot(lam0,betaPCM2Mat(1:p-3),'k*','MarkerSize',5)

figure;
plot(outPCM2path.lamPath(2:25)/n,mean(sigmaPCM2Matpath(1:p-3,2:25)/2),'LineWidth',3)
grid on
box on
hold on
plot(outPCM1path.lamPath(2:25)/n,mean(sigmaPCM1Matpath(1:p-3,2:25)),'LineWidth',3)
xlabel('\lambda path')
ylabel('\sigma_i')
set(gca,'FontSize',20)
legend('LS','Huber','Location','NorthWest')


% Plot lam0 solutions for comparison
% Two selected sets
selSet1 = find(abs(betaPCM1Mat)>epsThresh);
selSet2 = find(abs(betaPCM2Mat)>epsThresh);

% Joint sets
jointSet = union(selSet1,selSet2);
setLen = length(jointSet);

% Figure to show the solution using the theoretical lambda
figure;
stem((1:setLen)-0.15,betaPCM2Mat(jointSet),'LineWidth',6);
hold on;
stem((1:setLen)+0.15,betaPCM1Mat(jointSet),'LineWidth',6);
legend('LS','Huber','Location','NorthWest')
grid on
ylabel('\beta_i')
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen))
set(gca,'XTickLabel',PredLabels(jointSet),'fontsize',18)
box on
xtickangle(90)
xlim([0 setLen+0.5])


figure;
stem((1:p)-0.5,betaPCM2Mat,'LineWidth',5);
hold on;
stem(1:p,betaPCM1Mat,'LineWidth',5);
legend('LS','Huber','Location','NorthWest')
grid on
xlabel('Taxon i')
ylabel('Selection frequency')
set(gca,'FontSize',20)
set(gca,'XTick',(1:p)-0.25)
set(gca,'XTickLabel',jointSet,'fontsize',18)
box on
xtickangle(90)
