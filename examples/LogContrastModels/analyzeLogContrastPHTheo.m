%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for soil pH/microbiome data
%
% Analyzes the regularization path for Huber and LS model and the 
% theoretical lambda solution
% 
% Reproduces Appendix D in Combettes & Mueller, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thresholding the numerical solutions
epsThresh = 1e-3;
betaPCM1Matpath(abs(betaPCM1Matpath(:))<epsThresh)=0;
betaPCM2Matpath(abs(betaPCM2Matpath(:))<epsThresh)=0;

figure;
plot(outPCM1path.lamPath/n,betaPCM1Matpath,'LineWidth',3)
grid on
box on
xlabel('\lambda path')
ylabel('\beta_i')
ylim([-0.6 0.6])
set(gca,'FontSize',20)
hold on
plot(lam0,betaPCM1Mat,'k*','MarkerSize',5)
title('Huber path')

figure;
plot(outPCM2path.lamPath/n,betaPCM2Matpath,'LineWidth',3)
grid on
box on
xlabel('\lambda path')
ylabel('\beta_i')
ylim([-0.6 0.6])
set(gca,'FontSize',20)
hold on
plot(lam0,betaPCM2Mat,'k*','MarkerSize',5)
title('LS path')

figure;
plot(outPCM2path.lamPath/n,mean(sigmaPCM2Matpath),'LineWidth',3)
grid on
box on
hold on
plot(outPCM1path.lamPath/n,mean(sigmaPCM1Matpath),'LineWidth',3)
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

% Taxa labels

% Select order names for the joint set
jointTaxa = tax_table(jointSet,5)
jointTaxaLabels = cell(setLen,1);
for i=1:setLen
    temp = strsplit(char(jointTaxa.order(i)),'__')
    jointTaxaLabels{i} = temp{2};
    if isempty(temp{2})
        jointTaxaLabels{i} = 'unknown';
    end
end

% Figure to show the selection probabilities of all selected taxa
figure;
stem((1:setLen)-0.15,betaPCM2Mat(jointSet),'LineWidth',6);
hold on;
stem((1:setLen)+0.15,betaPCM1Mat(jointSet),'LineWidth',6);
legend('LS','Huber','Location','NorthWest')
grid on
ylabel('\beta_i')
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen))
set(gca,'XTickLabel',jointTaxaLabels,'fontsize',18)
box on
xtickangle(90)
xlim([0.25 setLen+0.5])
