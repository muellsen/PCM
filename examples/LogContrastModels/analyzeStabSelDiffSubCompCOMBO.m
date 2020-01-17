%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analyze results from log-contrast modeling on COMBO data
% microbiome + covariates with general subcompositional constraints
% and standard log-contrast constraints
%
% Analyze the influence of the subcompositional constraint on stability
% selection of the variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stability profile with subcompositional constraints
selSet1 =find(selFreq1~=0);
selSet2 =find(selFreq2~=0);

jointSetC = union(selSet1,selSet2);

% Stability profile with standard log-contrast constraints
selSet1ones =find(selFreq1ones~=0);
selSet2ones =find(selFreq2ones~=0);

jointSetones = union(selSet1ones,selSet2ones);

jointSet = union(jointSetones,jointSetC); 

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
ylim([0 1])

figure;stem((1:setLen)-0.5,selFreq2ones(jointSet),'LineWidth',5);hold on;stem(1:setLen,selFreq1ones(jointSet),'LineWidth',5);
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
ylim([0 1])

% Difference in frequencies
figure;stem((1:setLen)-0.5,selFreq2(jointSet)-selFreq2ones(jointSet),'LineWidth',5);hold on;stem(1:setLen,selFreq1(jointSet)-selFreq1ones(jointSet),'LineWidth',5);
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
ylim([-0.5 0.5])








