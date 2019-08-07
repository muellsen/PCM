% Analyze the experimental design
power1Vec = [2,3/2];
nPow1 = length(power1Vec);

objFunCell = {'Lq','Huber'}
nObj = length(objFunCell);

penFunCell = {'','Lq','BerHu'}
nPen = length(penFunCell);

power2Vec = [2];
nPow2 = length(power2Vec);

regCell = {'','L1','L1s'}
nReg = length(regCell);

modes = {'homosced','hetsced'}
nModes = length(modes);

nExperiments = nObj*nPow1*nPen*nPow2*nReg*nModes;

cnt = 1;

minPredErrVec = zeros(nExperiments,1);
minAbsErrVec = zeros(nExperiments,1);

minEstErrVec = zeros(nExperiments,1);
minSuppErrVec = zeros(nExperiments,1);

minPredIndVec = zeros(nExperiments,1);
minEstIndVec = zeros(nExperiments,1);
minSuppIndVec = zeros(nExperiments,1);

titleCell = cell(nExperiments,1);

for i1 = 1:nPow1
    for i2=1:nObj
        for i3 = 1:nPen
            for i4 = 1:nPow2
                for i5 = 1:nReg
                    for i6 = 1:nModes
                        
                        if strcmp(penFunCell{i3},'') && strcmp(regCell{i5},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1))];
                            
                        elseif strcmp(penFunCell{i3},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Reg: ',regCell{i5},', Mode:', modes{i6}];
                            
                        elseif strcmp(regCell{i5},'')
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Pen: ' penFunCell{i3}, ' with q=',num2str(power2Vec(i4)), ', Mode:', modes{i6}];
                        else
                            titleString = ['Fit: ', objFunCell{i2} , ' with q=',num2str(power1Vec(i1)), ...
                                ', Pen: ' penFunCell{i3}, ' with q=',num2str(power2Vec(i4)), ', Reg: ',regCell{i5},', Mode:', modes{i6}];
                        end
                        
                        currBeta = betaCell{cnt};
                        currSigma = sigmaCell{cnt};
                        currTau = tauCell{cnt};
                        currOut = outCell{cnt};
                        currRuntime = runTimeVec(cnt);
                        %
                        %                                                 figure(cnt);
                        %                                                 semilogx(currOut.lamPath,currBeta,'LineWidth',1)
                        %                                                 hold on
                        %                                                 semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        %                                                 xlabel('Regularization path \lambda')
                        %                                                 ylabel('\beta_i')
                        %                                                 semilogx(currOut.lamPath,repmat([-1;1],1,length(currOut.lamPath)),'k--','LineWidth',1)
                        %                                                 xlim([3 20])
                        %                                                 ylim([-2.5 2.5])
                        %                                                 set(gca,'FontSize',15)
                        %                                                 grid on
                        %                                                 title(titleString)
                        %                                                 drawnow;
                        %                                                 saveas(gcf,['FigFocus2_',num2str(cnt)],'png')
                        %                                                 close gcf
                        %
                        %                         figure(cnt+100);
                        %                         semilogx(currOut.lamPath,currTau,'LineWidth',1)
                        %                         hold on
                        %                         semilogx(currOut.lamPath,currTau(1:nnzs,:),'LineWidth',5)
                        %                         xlabel('Regularization path \lambda')
                        %                         ylabel('\tau_i')
                        %                         set(gca,'FontSize',15)
                        %                         grid on
                        %                         title(titleString)
                        %                         drawnow;
                        %                         %saveas(gcf,['FigTau_',num2str(cnt)],'png')
                        
                        %                         figure(cnt+200);
                        %                         semilogx(currOut.lamPath,currSigma([1,1+n_1,1+n_1+n_2],:),'LineWidth',5)
                        %                         hold on
                        %                         if strcmp(modes{i6},'hetsced')
                        %                             legend('group 1','group 2','group 3','Location','NorthWest')
                        %                         end
                        %                         semilogx(currOut.lamPath,repmat([s1;s2;s3],1,length(currOut.lamPath)),'k--','LineWidth',5)
                        %                         xlabel('Regularization path \lambda')
                        %                         ylabel('\sigma_i')
                        %                         ylim([0 7])
                        %
                        %                         set(gca,'FontSize',15)
                        %                         grid on
                        %                         title(titleString)
                        %                         drawnow;
                        %                         saveas(gcf,['FigSigma_',num2str(cnt)],'png')
                        
                        % Compute oracle minima of prediction, estimation,
                        % and support
                        [minErr,minInd] = min(sum(abs(currBeta-repmat(betaTrue,1,length(currOut.lamPath))).^2)./n);
                        
                        minEstErrVec(cnt) = minErr;
                        minEstIndVec(cnt) = minInd;
                        
                        [minErr,minInd] = min(sum(abs(X*currBeta-X*repmat(betaTrue,1,length(currOut.lamPath))).^2)./n);
                        minPredErrVec(cnt) = minErr;
                        minPredIndVec(cnt) = minInd;
                        
                        minAbsErrVec(cnt) = min(sum(abs(X*currBeta-X*repmat(betaTrue,1,length(currOut.lamPath))))./n);
                        
                        
                        titleCell{cnt} = titleString;
                        
                        cnt = cnt+1;
                    end
                end
            end
        end
    end
end

% Top best fits

[sortedVals,sortedInds] = sort(minEstErrVec,'ascend');

titleCell(sortedInds)
bestInd = sortedInds(2);

currBetaMat = betaCell{bestInd};
currBeta1 = currBetaMat(:,minPredIndVec(bestInd));
hvarInds = 1:n_1;

figure;
plot(X*currBeta1,Y,'.','MarkerSize',30)
grid on
hold on
plot(X(hvarInds,:)*currBeta1,Y(hvarInds),'r.','MarkerSize',30)
xlabel('Predicted response')
ylabel('True response')
set(gca,'FontSize',20);
legend({'Low variance','High variance'},'FontSize',20)
title(titleCell{bestInd})
xlim([-12 12])
ylim([-12 12])
saveas(gcf,['FigBestPred_',num2str(bestInd)],'png')

%
bestInd = sortedInds(9);

currBetaMat = betaCell{bestInd};
currBeta2 = currBetaMat(:,minPredIndVec(bestInd));

figure;
plot(X*currBeta2,Y,'.','MarkerSize',30)
grid on
hold on
plot(X(hvarInds,:)*currBeta2,Y(hvarInds),'r.','MarkerSize',30)
xlabel('Predicted response')
ylabel('True response')
set(gca,'FontSize',20);
legend({'Low variance','High variance'},'FontSize',20)
title(titleCell{bestInd})
xlim([-12 12])
ylim([-12 12])
saveas(gcf,['FigBestPred_',num2str(bestInd)],'png')

% Oracle
figure;
plot(X*betaTrue,Y1,'.','MarkerSize',30)
grid on
hold on
plot(X(hvarInds,:)*betaTrue,Y(hvarInds),'r.','MarkerSize',30)
xlabel('Predicted response')
ylabel('True response')
set(gca,'FontSize',20);
legend({'Low variance','High variance'},'FontSize',20)
title('Oracle predictor')
xlim([-12 12])
ylim([-12 12])
saveas(gcf,'FigBestOracle','png')







