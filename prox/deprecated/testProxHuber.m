% Test proximity operator of Huber function
close all

% Rho parameter
rho = 1.0;
gamma = 1;

% Create test points
p = 2;
n = 1e3;
const = 10;

testPoints = -const+2*const.*rand(p,n);

proxPoints = zeros(p,n);
caseMat = zeros(p-1,n);

% Five colors for five cases
colorVec={'b','g','m','y','k','r'};

for i=1:n
    
    eta = testPoints(1,i);
    y = testPoints(2:end,i);
    
    [etaProx,yProx,caseInds] = proxgpHuber(eta,y,rho,gamma);
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
        set(gca,'FontSize',30)
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
            xlim([-const const])
            ylim([-const const])
            set(gca,'FontSize',30)
            
            %drawnow
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
            xlim([-const const])
            ylim([-const const])
            zlim([-const const])
            
            figure(25)
            plot3([eta,eta-etaProx],[y(1),y(1)-yProx(1)],[y(2),y(2)-yProx(2)],[colorVec{caseInd},'-'],'LineWidth',5)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-const const])
            ylim([-const const])
            zlim([-const const])
            set(gca,'FontSize',30)
            
            figure(27)
            plot3(eta-etaProx,y(1)-yProx(1),y(2)-yProx(2),[colorVec{caseInd},'.'],'MarkerSize',20)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-100 100])
            ylim([-1 1])
            zlim([-1 1])
            set(gca,'FontSize',30)
            
            drawnow
        end
    end
    
    
end
