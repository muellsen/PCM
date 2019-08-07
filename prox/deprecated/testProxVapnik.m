% Test proximity operator of Vapnik function
close all

% Rho parameter
eps = 1;
gamma = 1;

% Create test points
p = 3;
n = 1e5;

testPoints = -2+4.*rand(p,n);

proxPoints = zeros(p,n);
caseMat = zeros(p-1,n);

% Five colors for five cases
colorVec={'b','g','m','y','k','g','m','r'};

for i=1:n
    
    eta = testPoints(1,i);
    y = testPoints(2:end,i);
    
    [etaProx,yProx,caseInds] = proxgVapnik(eta,y,eps,gamma);
    caseMat(:,i) = caseInds;
    proxPoints(1,i) = etaProx;
    proxPoints(2:end,i) = yProx;
    
    caseInd = caseMat(1,i);
    
    if p==2
        
        figure(23)
        plot(eta,y,[colorVec{caseInd},'.'],'MarkerSize',50)
        hold on
        grid on
        box on
        xlabel('eta')
        ylabel('y')
        drawnow
        %
        %     figure(24)
        %     plot(eta-etaProx,y-yProx,[colorVec{caseInd},'.'],'MarkerSize',50)
        %     hold on
        %     grid on
        %     box on
        %     xlabel('eta')
        %     ylabel('y')
        %     drawnow
        
        if caseInd>0
            figure(25)
            plot([eta,eta-etaProx],[y,y-yProx],[colorVec{caseInd},'-'],'LineWidth',5)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y')
            xlim([-2 2])
            ylim([-2 2])
            drawnow
        end
    elseif p==3
        
        if caseInd>0
            
            figure(23)
            plot3(eta,y(1),y(2),[colorVec{caseInd},'.'],'MarkerSize',50)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            xlim([-2 2])
            ylim([-2 2])
            zlim([-2 2])
            
            
            %drawnow
            %
            %     figure(24)
            %     plot(eta-etaProx,y-yProx,[colorVec{caseInd},'.'],'MarkerSize',50)
            %     hold on
            %     grid on
            %     box on
            %     xlabel('eta')
            %     ylabel('y')
            %     drawnow
            
            figure(25)
            plot3([eta,eta-etaProx],[y(1),y(1)-yProx(1)],[y(2),y(2)-yProx(2)],[colorVec{caseInd},'-'],'LineWidth',5)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-2 2])
            ylim([-2 2])
            zlim([-2 2])
            set(gca,'FontSize',30)
            
            figure(27)
            plot3(eta-etaProx,y(1)-yProx(1),y(2)-yProx(2),[colorVec{caseInd},'.'],'MarkerSize',20)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-2 2])
            ylim([-2 2])
            zlim([-2 2])
            set(gca,'FontSize',30)
            
            drawnow
        end
    end
end
