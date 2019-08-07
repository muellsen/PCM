% Make nice plots for Patrick

figure(1)
hold off
plot(lassoInfo.Lambda(:,end:-1:1),betaLassoMat')
title('Lasso','FontSize',20)
drawnow
hold on
plot(lassoInfo.Lambda(:,end:-1:1),betaLassoMat(1:10,:)','LineWidth',5)
xlabel('Regularization parameter \lambda')
ylabel('\beta_j estimates')
set(gca,'FontSize',20)
grid on
box on

pathFac = (8*sqrt(sum(Y.^2))*sqrt(n-1))/(n);

figure(5)
hold off
plot(outPCM.lamPath/pathFac,betaPCMMat')
title('Robust concomitant estimation','FontSize',20)
drawnow
hold on
plot(outPCM.lamPath/pathFac,betaPCMMat(1:nnzs,:)','LineWidth',5)
xlabel('Regularization parameter \lambda')
ylabel('\beta_j estimates')
set(gca,'FontSize',20)
grid on
box on

trueScales = repmat([sig1(1);sig2(1);sig3(1)],1,length(outPCM.lamPath))

figure(51)
hold off
plot(outPCM.lamPath/pathFac,sigmaPCMMat','LineWidth',5)
title('Robust concomitant estimation','FontSize',20)
drawnow
hold on
plot(outPCM.lamPath/pathFac,trueScales,'k--','LineWidth',2)
xlabel('Regularization parameter \lambda')
ylabel('Concomitant scales \sigma_i')
set(gca,'FontSize',20)
grid on
box on

res1 = (X(1:n_1,:)*betaPCMMat-repmat(Y(1:n_1),1,length(outPCM.lamPath)));
res2 = (X(n_1+1:n_1+n_2,:)*betaPCMMat-repmat(Y(n_1+1:n_1+n_2),1,length(outPCM.lamPath)));
res3 = (X(n_1+n_2+1:end,:)*betaPCMMat-repmat(Y(n_1+n_2+1:end),1,length(outPCM.lamPath)));

stdEst1b = sqrt(sum(res1.^2)/n_1)
stdEst2b = sqrt(sum(res2.^2)/n_2)
stdEst3b = sqrt(sum(res3.^2)/n_3)

figure(501)
hold off
plot(outPCM.lamPath/pathFac,[stdEst1b;stdEst2b;stdEst3b],'LineWidth',5)
title('Robust concomitant estimation','FontSize',20)
drawnow
hold on
plot(outPCM.lamPath/pathFac,trueScales,'k--','LineWidth',2)
xlabel('Regularization parameter \lambda')
ylabel('Concomitant scales \sigma_i')
set(gca,'FontSize',20)
grid on
box on

% Mean shift refit and reduced sigma estimation

% Estimate sigma by discarding outliers

currY = Y(1:n_1); 
currX = X(1:n_1,:);
lenPath = length(outPCM.lamPath);
gamMat1 = zeros(n_1,lenPath);
stdEst1 = zeros(1,lenPath);

for i=1:lenPath
    i
    gamMat1(:,i) = MeanShiftrefit(currY,currX,betaPCMMat(:,i),sigmaPCMMat(1,i),proxopts.rho1,1);
    currInlInds = abs(gamMat1(:,i))<1e-1;
    currRes = (currX(currInlInds,:)*betaPCMMat(:,i)-currY(currInlInds));
    stdEst1(i) = sqrt(sum(currRes.^2)/(sum(currInlInds)))
end

currY = Y(n_1+1:n_1+n_2); 
currX = X(n_1+1:n_1+n_2,:);
gamMat2 = zeros(n_2,lenPath);
stdEst2 = zeros(1,lenPath);

for i=1:lenPath
    i
    gamMat2(:,i) = MeanShiftrefit(currY,currX,betaPCMMat(:,i),sigmaPCMMat(2,i),proxopts.rho1,1);
    currInlInds = abs(gamMat2(:,i))<1e-1;
    currRes = (currX(currInlInds,:)*betaPCMMat(:,i)-currY(currInlInds));
    stdEst2(i) = sqrt(sum(currRes.^2)/(sum(currInlInds)))
end

currY = Y(n_1+n_2+1:end); 
currX = X(n_1+n_2+1:end,:);
gamMat3 = zeros(n_3,lenPath);
stdEst3 = zeros(1,lenPath);

for i=1:lenPath
    i
    gamMat3(:,i) = MeanShiftrefit(currY,currX,betaPCMMat(:,i),sigmaPCMMat(3,i),proxopts.rho1,1);
    currInlInds = abs(gamMat3(:,i))<1e-1; 
    currRes = (currX(currInlInds,:)*betaPCMMat(:,i)-currY(currInlInds));
    stdEst3(i) = sqrt(sum(currRes.^2)/(sum(currInlInds)))
end


figure(501)
hold off
plot(outPCM.lamPath/pathFac,[stdEst1;stdEst2;stdEst3],'LineWidth',5)
title('Robust concomitant estimation','FontSize',20)
drawnow
hold on
plot(outPCM.lamPath/pathFac,trueScales,'k--','LineWidth',2)
xlabel('Regularization parameter \lambda')
ylabel('Concomitant scales \sigma_i')
set(gca,'FontSize',20)
grid on
box on


% Oracle estimation
stdEst1o = sqrt((sum(((X(inlInds1,:)*repmat(betaTrue,1,length(outPCM.lamPath))-repmat(Y(inlInds1),1,length(outPCM.lamPath)))).^2))/(length(inlInds1)))
stdEst2o = sqrt((sum(((X(inlInds2,:)*repmat(betaTrue,1,length(outPCM.lamPath))-repmat(Y(inlInds2),1,length(outPCM.lamPath)))).^2))/(length(inlInds2)))
stdEst3o = sqrt((sum(((X(inlInds3,:)*repmat(betaTrue,1,length(outPCM.lamPath))-repmat(Y(inlInds3),1,length(outPCM.lamPath)))).^2))/(length(inlInds3)))

figure(501)
hold off
plot(outPCM.lamPath/pathFac,[stdEst1o;stdEst2o;stdEst3o],'LineWidth',5)
title('Robust concomitant estimation','FontSize',20)
drawnow
hold on
plot(outPCM.lamPath/pathFac,trueScales,'k--','LineWidth',2)
xlabel('Regularization parameter \lambda')
ylabel('Concomitant scales \sigma_i')
set(gca,'FontSize',20)
grid on
box on




















