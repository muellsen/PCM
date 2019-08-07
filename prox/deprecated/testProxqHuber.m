% Test proximity operator of Huber function
close all

% Parameter
rho = 1.0;
gamma = 1;

qPower = 2;
alpha = 1;
kappa = 1;

% Create test points
p = 3;
n = 1e3;
const = 10;

testPoints = -const+2*const.*rand(p,n);

proxPoints = zeros(p,n);
caseMat = zeros(1,n);

% Five colors for five cases
colorVec={'b','g','m','y','k','r'};

for i=1:n
    
    eta = testPoints(1,i);
    y = testPoints(2:end,i);
    
    %[etaProx,yProx,caseInds] = proxgqHuber(eta,y,rho,alpha,qPower,gamma);
    %[etaProx,yProx,caseInds] = proxgqBerHu(eta,y,rho,alpha,kappa,qPower,gamma);
    %[etaProx2,yProx2,caseInds2] = proxgpBerHu(eta,y,rho,gamma);
    tic
    %[etaProx,yProx,caseInds] = proxgqLq(eta,y,kappa,alpha,qPower,gamma);
    [etaProx,yProx,caseInds] = proxgqHuber(eta,y,rho,alpha,qPower,gamma);
    time1=toc
    tic
    %[etaProx2,yProx2,caseInds2] = proxg2Lq(eta,y,kappa,alpha,gamma);
    [etaProx2,yProx2,caseInds2] = proxg2Huber(eta,y,rho,alpha,gamma);
    
    time2=toc
    
    etaProx2-etaProx
    yProx2-yProx
    if abs(sum(yProx2-yProx))>1e-6
        error('Prox diff')
    end
    
    caseMat(1,i) = caseInds;
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
        set(gca,'FontSize',30)
        
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
        
        xlim([-const const])
        ylim([-const const])
        zlim([-const const])
        set(gca,'FontSize',30)
        
        drawnow
    end
    
end
