
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log-contrast modeling for soil pH/microbiome data
%
% Analyze stability selection + refitting for Huber and LS model
%
% Reproduces (up to variation in stability selection) Section 4.2 in
% Combettes & Mueller, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Two selected sets
selSet1 = find(selFreq1>0.1);
selSet2 = find(selFreq2>0.1);

% Joint sets
jointSet = union(selSet1,selSet2);
setLen = length(jointSet);

% Idenitfy taxonomic names at lowest available rank
jointTaxa = correctTaxTable(jointSet,:);
nTab = size(jointTaxa,2);

jointTaxaLabels = cell(setLen,1);
for i=1:setLen
    for j=nTab-1:-1:2
        currEntry = table2array(jointTaxa(i,j));
        temp = strsplit(char(currEntry),'__')
        if ~isempty(temp{2})
            jointTaxaLabels{i} = temp{2};
            break;
        end
    end
end

% Figure to show the selection probabilities of all selected taxa
figure;stem((1:setLen)-0.25,selFreq2(jointSet),'LineWidth',5);
hold on;stem((1:setLen)+0.25,selFreq1(jointSet),'LineWidth',5);
legend('LS','Huber','Location','NorthWest')
grid on
xlabel('Taxon i')
ylabel('Selection frequency')
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen))
set(gca,'XTickLabel',jointTaxaLabels,'fontsize',18)
box on
xtickangle(90)
xlim([0.25 setLen+0.5])

% Threshold
th_stab = 0.7;

selTop1 = find(selFreq1>=th_stab)
selTop2 = find(selFreq2>=th_stab)

% Reoptimize over the subset
X_sel1 = X(:,selTop1);

[n,p] = size(X_sel1);

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

pcmopts.abstol = 1e-7;
pcmopts.lamPath = 0;% no penalty;
pcmopts.gamma = 1;

t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(X_sel1, Y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% Reoptimize over the subset for LS model
X_sel2 = X(:,selTop2);
[n,p] = size(X_sel2);

% Linear constraint for log-contrast model
Ceq = ones(1,p); % Standard log-constrast constraint
rhsvec = zeros(size(Ceq,1),1);

% Optimization model
% Algorithmic parameters
clear pcmopts;

% Optimization model
power1 = 2;
objFun = 'Lq';
penFun = 'L1';
regFun = 'Ceq';

pcmopts.qPower1 = power1;
pcmopts.objFun = objFun;
pcmopts.fitLin = 1/2;
pcmopts.penFun = penFun;
pcmopts.regFun = regFun;
pcmopts.Ceq = Ceq;
pcmopts.rhsvec = rhsvec;

pcmopts.abstol = 1e-7;
pcmopts.lamPath = 0;% no penalty
pcmopts.gamma = 1;

t1=now;
[betaPCM2Mat, sigmaPCM2Mat,funPCM2Mat,outPCM2] = pcmC2(X_sel2, Y_cent, pcmopts);
t2=now;
timePCM2 = (t2-t1)*(60*60*24)

% Identify union of identified variables
jointSelTop = union(selTop1,selTop2);

% Identify OTU index for both selected sets in joint set
[~,indTop1] = intersect(jointSelTop,selTop1);
[~,indTop2] = intersect(jointSelTop,selTop2);

% Idenitfy taxonomic names at lowest available rank
topTaxa = correctTaxTable(jointSelTop,:);

setLen=size(topTaxa,1);

topLabels = cell(setLen,1);
for i=1:setLen
    for j=nTab-1:-1:2
        currEntry = table2array(topTaxa(i,j));
        temp = strsplit(char(currEntry),'__')
        if ~isempty(temp{2})
            topLabels{i} = temp{2};
            break;
        end
    end
end

% Insert by hand (deprecated)
% topNames = {'Ellin6513-1','RB41','Balneimonas-1','Rubrobacter','Ellin6513-2','Koribacteraceae','Balneimonas-2'}

% Figure to show the beta values for the selected variables
figure;
stem((indTop2)-0.15,betaPCM2Mat,'LineWidth',10);
hold on;
stem((indTop1)+0.15,betaPCM1Mat,'LineWidth',10);
legend('LS','Huber','Location','NorthWest')
grid on
ylabel('\beta_i','FontSize',40)
set(gca,'FontSize',20)
set(gca,'XTick',(1:setLen))
set(gca,'XTickLabel',topLabels,'fontsize',18)
box on
xtickangle(90)
xlim([0.5 setLen+0.5])
ylim([-0.4 0.4])

% Analyze Huber model refit
currOut = outPCM1;
currBeta = betaPCM1Mat;
currSigma = sigmaPCM1Mat;

currInd = 1;

res = (Y_cent-X_sel1*currBeta(:,currInd));
eta = currSigma(:,currInd);
[fun,quadInds,linInds] = HuberFun(eta,res,currOut.opts.rho1,currOut.opts.fitLin(1),currOut.opts.qPower1);

% Figure to show the selection probabilities of all selected taxa
figure;
plot(y_bar+X_sel1(quadInds,:)*betaPCM1Mat,Y(quadInds),'.','MarkerSize',30);
hold on;
plot(y_bar+X_sel1(linInds,:)*betaPCM1Mat,Y(linInds),'*','MarkerSize',15);
legend('Inliers','Outliers','Location','NorthWest')
grid on
xlabel('pH predictions','FontSize',20)
ylabel('pH measurements','FontSize',40)
set(gca,'FontSize',20)
h=line([3,9],[3,9],'Color','black','LineStyle','--','Linewidth',2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Summary of the model fit for the Huber model

% Mean absolute deviation
MAD_all = mean(abs(res));
MAD_in = mean(abs(res(quadInds)));
MAD_out = mean(abs(res(linInds)));

% LS deviation
LS_all = sqrt(sum(res.^2)/length(res));
LS_in = sqrt(sum((res(quadInds).^2)/length(quadInds)));
LS_out = sqrt(sum((res(linInds).^2))/length(linInds));

% R^2
R2 = (corr(Y_cent,X_sel1*currBeta(:,currInd))).^2

% Analyze LS model refit

currOut = outPCM2;
currBeta = betaPCM2Mat;
currSigma = sigmaPCM2Mat;

currInd = 1;
res2 = (Y_cent-X_sel2*currBeta(:,currInd));

% Summary of the model fit for the LS model

% Mean absolute deviation
MAD_all2 = mean(abs(res2));

% LS deviation
LS_all2 = sqrt(sum(res2.^2)/length(res2));

% R^2
R2_2 = (corr(Y_cent,X_sel2*currBeta(:,currInd))).^2



