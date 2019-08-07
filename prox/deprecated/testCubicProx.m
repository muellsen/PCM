
% Test proximity operator of cubic function
%close all

% Rho parameter
eps = 1;
gamma = 1;

% Create test points
p = 3;
n = 1e3;

const = 100;
testPoints = -const+2*const.*rand(p,n);

proxPoints = zeros(p,n);
caseMat = zeros(p-1,n);

% Five colors for five cases
colorVec={'b','g','k','y','m','r','w','r'};

for i=1:n
    
    eta = testPoints(1,i);
    y = testPoints(2:end,i);
    
    % Set all parameters
    v = zeros(p-1,1);
    Yvec = zeros(p-1,1);
    
    [etaProx,yProx] = proxgTREX(eta,y,v,Yvec)
    eta_y = [eta;y];
    alpha = 1;
    gamma = 1;
    [temp,caseInds] = proxgTREX2(eta_y,v,Yvec,alpha,gamma);
    etaProx2 = temp(1);
    yProx2 = temp(2:end);
    
    
    %caseInds = 1;
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
            xlim([-const const])
            ylim([-const const])
            drawnow
        end
    elseif p==3
        
        if caseInd>0
            
            figure(50)
            plot3(eta,y(1),y(2),[colorVec{caseInd},'.'],'MarkerSize',50)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            xlim([-const 0])
            ylim([-const/4 const/4])
            zlim([-const/4 const/4])
                        set(gca,'FontSize',30)

            
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
            
            plot3([eta,eta-etaProx],[y(1),y(1)-yProx(1)],[y(2),y(2)-yProx(2)],[colorVec{caseInd},'-'],'LineWidth',5)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-const 0])
            ylim([-const/4 const/4])
            zlim([-const/4 const/4])
            
            set(gca,'FontSize',30)
            
            figure(42)
            plot3(eta-etaProx,y(1)-yProx(1),y(2)-yProx(2),[colorVec{caseInd},'.'],'MarkerSize',20)
            hold on
            grid on
            box on
            xlabel('eta')
            ylabel('y_1')
            zlabel('y_2')
            
            xlim([-const 0])
            ylim([-const/4 const/4])
            zlim([-const/4 const/4])
            
            set(gca,'FontSize',30)
            
            drawnow
        end
    end
end
