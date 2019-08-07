% Test proximity operator of Huber function
close all

% Rho parameter
rho = 0.5;
gamma = 1;

% Create test points
p = 100;
n = 1e4;

testPoints = -2+4.*randn(p,n);

proxPoints = zeros(p,n);
caseMat = zeros(p-1,n);

for i=1:n
    
    eta = testPoints(1,i);
    y = testPoints(2:end,i);
    
    [etaProx,yProx,caseInds] = proxgpHuber(eta,y,rho,gamma);
    caseMat(:,i) = caseInds;
    proxPoints(1,i) = etaProx;
    proxPoints(2:end,i) = yProx;
    
end

dPoints = (testPoints(:,2:end)-testPoints(:,1:end-1));
dProxs = (proxPoints(:,2:end)-proxPoints(:,1:end-1));
dProxPoints = ((testPoints(:,1:end-1)-proxPoints(:,1:end-1))...
              -(testPoints(:,2:end)-proxPoints(:,2:end)));
                  
          
normP = sum(dPoints.^2);
normDiffProx = sum((dProxs).^2);    
normDiffPProx = sum(dProxs.*dPoints);
    
proxDiffVec = (normDiffPProx-normDiffProx)
    
figure;
plot(proxDiffVec,'.','MarkerSize',50)
grid on
box on
set(gca,'FontSize',30)








