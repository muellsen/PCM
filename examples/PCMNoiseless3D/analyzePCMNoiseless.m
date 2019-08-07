% Analyze PCM solver for low dimensional example loss function on the same data set with
% extension to outlier removal

load NoiseFreeOutlier.mat

cnt = 1;

for i1 = 1:nPow1
    for i2=1:nObj
        for i3 = 1:nLB
            for i4 = 1:nData
                for i5 = 1:nReg
                    for i6 = 1:nModes
                        
                        % With or without outlier
                        if i4==2
                            outlString = 'Yes';
                        else
                            outlString = 'No';
                        end
                        
                        
                        titleString = ['Fit: ', objFunCell{i2},', Reg: ',regCell{i5},', Outlier: ',outlString,', Mode: ', modes{i6},', LB: ', num2str(lbVec(i3))];
                        
                        currBeta = betaCell{cnt};
                        currSigma = sigmaCell{cnt};
                        currTau = tauCell{cnt};
                        currOut = outCell{cnt};
                        currRuntime = runTimeVec(cnt);
                        
                        figure(cnt);
                        plot(currOut.lamPath,currBeta,'LineWidth',5)
                        hold on
                        %semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        xlabel('Regularization parameter \alpha_1','FontSize',30)
                        ylabel('\beta_i','FontSize',30)
                        grid on
                        %title(titleString)
                        plot(currOut.lamPath,repmat([-const;const],1,length(currOut.lamPath)),'k--','LineWidth',5)
                        xlim([8e-2 10])
                        set(gca,'FontSize',20)
                        
                        drawnow;
                        saveas(gcf,['Fig3Dbeta_',num2str(cnt)],'png')
                        
                        figure(cnt+23);
                        plot(currOut.lamPath,currSigma,'LineWidth',5)
                        hold on
                        %semilogx(currOut.lamPath,currBeta(1:nnzs,:),'LineWidth',5)
                        xlabel('Regularization parameter \alpha_1','FontSize',30)
                        ylabel('\sigma_i','FontSize',30)
                        plot(currOut.lamPath,repmat([3;0],1,length(currOut.lamPath)),'k--','LineWidth',2)
                        
                        grid on
                        %title(titleString)
                        xlim([8e-2 10])
                        ylim([-0.5 3.5])
                        set(gca,'FontSize',20)
                        saveas(gcf,['Fig3Dsigma_',num2str(cnt)],'png')
                        %saveas(gcf,['Fig3DsigmaZoom2_',num2str(cnt)],'png')
                        drawnow;
                        
                        cnt = cnt+1;
                        
                    end
                end
            end
        end
    end
end





