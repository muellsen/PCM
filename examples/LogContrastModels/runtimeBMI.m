% Comparison of the runtime of exact proximal scheme and the coordinate
% descent (with line search) method (Shi et al. 2016)
% for sparse log-contrast regression with joint scale estimation
% on the BMI example without covariates. 
% (Appendix C in Combettes & Müller, 2019)

% Ensure that you are in the LogContrastModels folder

% Load BMI data in MATLAB format
load('data/BMIData_subset.mat')

% Input:
% x: n by p design matrix
p = nFiltSpecies;
n = nSamples;

% Total sum normalization + 0.5 pseudocount
x_c = (FiltSpeciesCntData+0.5)./repmat(sum((FiltSpeciesCntData+0.5)),nFiltSpecies,1);
x = log(x_c)';

% y: response of length n
y = BMIfilt;

% Center X and Y
y_bar = mean(y);
y_cent = y-y_bar;

loc_x = mean(x);
sca_x = std(x);
b2 = 1./sca_x;
diagB = diag(b2);
x_cent = (x - ones(n,1)*loc_x)*diagB;

%x2 = clr(x_c,0)';
%x_cent = x2;

x3 = clr(FiltSpeciesCntData,1)';

% constr: linear constraint, for compositional data, 
% here groups of four phyla (CFiltered is a px4 matrix)
constr = CFiltered;

% Linear constraint for log-contrast model
Ceq = constr'; % sub-compositional coherence
rhsvec = zeros(size(Ceq,1),1);

% Solution tolerance
absTol = 1e-8;

%% Perspective log-contrast model (Combettes and Mueller, 2019)

% Algorithmic parameters
clear pcmopts;

% Set-up experimental design
power1 = 2;
objFun = 'Lq'; % (LS for objective function) 
penFun = 'L1'; % (L1 for constraint)
regFun = 'Ceq'; % Linear constraint

% Specific model setup; here LS + L1 + simple linear constraint
pcmopts.qPower1 = power1;
pcmopts.objFun = objFun;
pcmopts.fitLin = 1/2;
pcmopts.penFun = penFun;
pcmopts.regFun = regFun;

% Constraint definition
pcmopts.Ceq = Ceq;
pcmopts.rhsvec = rhsvec;

% Algorithmic parameters
pcmopts.abstol = absTol;
pcmopts.gamma = 0.5;

% Regularization path length
pcmopts.lenPath = 50;

% Solve the model using the DR method 
t1=now;
[betaPCMMat, sigmaPCMMat,funPCMMat,outPCM] = pcmC2(x_cent, y_cent, pcmopts);
t2=now;
timePCM = (t2-t1)*(60*60*24)

%% Scaled constraint Lasso (with line search on scale) (Shi et al. 2016)

% Algorithmic parameters
clear concomopts;

concomopts.Ceq = Ceq;
concomopts.abstol = absTol;
concomopts.tol_sig = 10*absTol;

% Path length
concomopts.lenPath = 50;

% Solve the model using the method of Shi et al. (2016) 
t1=now;
[betaConLMat, sigmaConLMat,outConL] = concomlasso(x_cent, y_cent, concomopts);
t2=now;
timeConcom = (t2-t1)*(60*60*24);


%% Plot solutions and differences between the two methods for solving 
%% the LS solution with joint scale estimation

figure;stem(betaPCMMat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',20);
xlabel('\beta_i')
grid on
box on
title('Solution across \lambda-path (Douglas-Rachford (this work))','FontSize',15)

figure;stem(betaConLMat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',20);
xlabel('\beta_i')
grid on
box on
title('Solution across \lambda-path (Shi et al 2016)','FontSize',15)

figure;stem(betaConLMat-betaPCMMat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',20);
xlabel('Difference in \beta_i')
grid on
box on
title('Difference in solutions across \lambda-path','FontSize',15)

figure;
plot(outPCM.lamPath/n,(betaConLMat-betaPCMMat),'LineWidth',3)
ax = gca;
xlabel('Regularization path \lambda')
ylabel('Differnence in \beta_i')
ylim([-2.2e-6 2.2e-6])
grid on
box on
set(ax,'FontSize',20);
title('Difference in solutions across \lambda-path','FontSize',15)

figure;
plot(outConL.lamPath,betaConLMat,'LineWidth',3)
ax = gca;
set(ax,'FontSize',20);
xlabel('Regularization path \lambda')
grid on
box on
title('Solution path (Shi et al 2016)','FontSize',15)

figure;
plot(outPCM.lamPath/n,betaPCMMat,'LineWidth',3)
ax = gca;
set(ax,'FontSize',20);
xlabel('Regularization path \lambda')
grid on
box on
title('Solution path (Douglas-Rachford)','FontSize',15)

figure;
plot(outPCM.lamPath/n,mean(sigmaPCMMat/2),'LineWidth',8)
hold on
plot(outConL.lamPath,mean(sigmaConLMat/2),'--','LineWidth',8)
ax = gca;
set(ax,'FontSize',20);
xlabel('Regularization path \lambda')
grid on
box on
title('Scale estimate','FontSize',15)
legend('Douglas-Rachford (this work)','Shi et al. (2016)','Location','NW')

figure;
plot(outPCM.lamPath/n,mean(sigmaPCMMat)-mean(sigmaConLMat),'LineWidth',4)
ax = gca;
set(ax,'FontSize',20);
xlabel('Regularization path \lambda')
ylim([-1e-6 1e-6])
grid on
box on
title('Difference in scale estimates','FontSize',15)

