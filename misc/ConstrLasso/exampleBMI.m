% Run the BMI example without covariates from Shi et al. 2016

% Load all data in MATLAB format
load allData.mat;

%% input:
% x: n by p design matrix

% Total sum normalization + 0.5 pseudocount
x_c = (FiltSpeciesCntData+0.5)./repmat(sum((FiltSpeciesCntData+0.5)),nFiltSpecies,1);
x = log(x_c)';

% y: response of length n
y = BMIfilt;

% constr: linear constraint, for compositional data, constr=ones(p,1)/sqrt(p); for no constraint, constr=zeros(p,1)
% level: level of confidence interval
% tol: tolerance level in the computation
% maxiter: the maximum number of iteration in solving the regularized estimator using coordinate descent method

constr = CFiltered;
level = 0.95;
tol = 1e-8;
maxiter = 200;

tol_sig = 1e-9;

%% compute the projection matrix according to the constraints
[constr2,S,V] = svd(constr);
constr = constr2(:,1:size(S,2));
Pc = constr*constr';
x = x-x*Pc;

%% get regularized estimator
method = 1; % using coordinate descent
%method = 2; % using CVX
[res.bet_n res.int res.sigma res.lam0]=sslr(y, x, constr, tol, tol_sig,maxiter, method);

%% get confidence interval
% for faster computation, we suggest you use CVX when p>2*n
[res.bet_u res.CI_l res.CI_u res.M]=cdmm_ci(y, x, constr, res.bet_n, res.int, res.sigma, res.lam0, level, tol, maxiter); % using coodinate descent
%[res.bet_u res.CI_l res.CI_u res.M]=cvx_ci(y, x, constr, res.bet_n, res.int, res.sigma, res.lam0, level, tol); % using CVX

%% output
% bet_n: regularized estimator
% int: intercept
% sigma: estimated noise level
% bet_u: debiased estimator
% [CI_l, CI_u]: confidence interval

%% Plot results

figure;stem(res.bet_n)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',14);
grid on
box on
title('Sparse solution')

figure;stem(res.bet_u)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',14);
grid on
box on
title('Debiased solution')

% Computing R^2 statistics for the sparse estimator
ss_res_n = sum((y-(mean(y)+x*res.bet_n)).^2)
ss_tot_n = sum((y-mean(y)).^2)

R2_n = 1 - ss_res_n/ss_tot_n

corBMI_n = corr(x*res.bet_n,y)

figure;plot(y,mean(y)+x*res.bet_n,'.','MarkerSize',30)
ax = gca;
set(ax,'FontSize',14);
grid on
box on
title(['Sparse R^2 = ',num2str(R2_n),' \rho = ',num2str(corBMI_n)])
xlabel('Observed BMI')
ylabel('Sparse Fitted BMI')

% Computing R^2 statistics for the sparse estimator
ss_res_u = sum((y-(mean(y)+x*res.bet_u)).^2)
ss_tot_u = sum((y-mean(y)).^2)

R2_u = 1 - ss_res_u/ss_tot_u

corBMI_u = corr(x*res.bet_u,y)

figure;plot(y,mean(y)+x*res.bet_u,'.','MarkerSize',30)
ax = gca;
set(ax,'FontSize',14);
grid on
box on
title(['Debiased R^2 = ',num2str(R2_u),' \rho = ',num2str(corBMI_u)])
xlabel('Observed BMI')
ylabel('Debiased Fitted BMI')


% Comparison with PCM

p = nFiltSpecies;
n = nSamples;

% Center X and Y
y_bar = mean(y);
y_cent = y-y_bar;

loc_x = mean(x);
sca_x = std(x);
b2 = 1./sca_x;
diagB = diag(b2);
x_c = (x - ones(n,1)*loc_x)*diagB;
x_cent = x_c;


% Theoretical lambda
options = optimset('Display','off');
kk = fsolve(@(k) (norminv(1-k/p))^4 + 2*((norminv(1-k/p))^2) - k, p/2, options);
lam0 = sqrt(2/n)*norminv(1-kk/p);


% Linear constraint for log-contrast model
Ceq = constr'; % sub-compositional coherence
rhsvec = zeros(size(Ceq,1),1);

% Algorithmic parameters
clear pcmopts;

% Set-up experimental design
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

pcmopts.abstol = 1e-5;
%pcmopts.lamPath = n*lam0;
pcmopts.gamma = 0.6;
pcmopts.lenPath = 30;

t1=now;
[betaPCM1Mat, sigmaPCM1Mat,funPCM1Mat,outPCM1] = pcmC2(x_cent, y_cent, pcmopts);
t2=now;
timePCM1 = (t2-t1)*(60*60*24)

% Algorithmic parameters
clear concomopts;
%concomopts.lamPath = lam0;
concomopts.Ceq = Ceq;
concomopts.abstol = 1e-5;
concomopts.lenPath = 30;

t1=now;
[betaConLMat, sigmaConLMat,outConL] = concomlasso(x_cent, y_cent, concomopts);
t2=now;
timeConcom1 = (t2-t1)*(60*60*24);

figure;stem(betaPCM1Mat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',14);
grid on
box on
title('Sparse solution (PCM)')

figure;stem(betaConLMat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',14);
grid on
box on
title('Sparse solution (Shi et al)')

figure;stem(betaConLMat-betaPCM1Mat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
ax.XTickLabel = generaCell(:,6);
ax.XTickLabelRotation = 90;
set(ax,'FontSize',14);
grid on
box on
title('Difference solution')

figure;
plot(mean(sigmaConLMat),betaConLMat)
ax = gca;
ax.XTick =  [1:nFiltSpecies];
set(ax,'FontSize',14);
grid on
box on
title('Solution path (Shi et al)')

figure;
plot(mean(sigmaPCM1Mat),betaPCM1Mat)
ax = gca;
set(ax,'FontSize',14);
grid on
box on
title('Solution path (PCM)')

figure;
plot(outPCM1.lamPath,mean(sigmaPCM1Mat),'LineWidth',3)
hold on
plot(outConL.lamPath*n,mean(sigmaConLMat),'LineWidth',3)
ax = gca;
set(ax,'FontSize',14);
grid on
box on
title('Scale estimate')
legend('PCM','Shi et al.')







