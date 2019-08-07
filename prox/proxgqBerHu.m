function [etaProx,yProx,caseInds] = proxgqBerHu(eta,y,rho,alpha,kappa,qPower,gamma)
% Perspective of generalized BerHu in (p+1) dimensions

% eta: concomitant
% y: input
% rho: BerHu's rho parameter
% alpha: constant in linear part
% kappa: Slope of the inner BerHu part
% qPower: qth power of the outer BerHu part
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

% Dual norm
q_s = qPower/(qPower-1);

% Delta(y/gamma)

ygamnorm = sqrt(sum((y/gamma).^2));
y_norm = sqrt(sum(y.^2));

Delta_ygam = max(ygamnorm-kappa,0) + 1/q_s*(max(ygamnorm-kappa,0))^q_s;

% Check four cases for prox calculation

% Case 1
if Delta_ygam<=(alpha-eta/gamma)/rho
    
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif y_norm>=gamma*kappa+rho*(eta-gamma*alpha) && Delta_ygam > (alpha-eta/gamma)/rho
    
    
    % Real part of the polynomial (for numerical reasons only)
    polR = @(t) real(rho*(eta-gamma*alpha+gamma*rho*(t-kappa+1/q_s.*(t-kappa).^q_s))*(1+(t-kappa).^(q_s-1))+gamma.*t - y_norm);
    
    t = fzero(polR,y_norm);
    %t2 = bisection(polR,0,10); % Alternative root finder
    %t-t2
    %pause
    fmin = polR(t);
    
    if fmin>1e-4
        error('Prox failed with 1e-4 tolerance')
    end
    
    % Least-square formulation for the root
    %polT = @(t) real((rho*(eta-gamma*alpha+gamma*rho*(t-kappa+1/q_s.*(t-kappa).^q_s))*(1+(t-kappa).^(q_s-1))+gamma.*t - y_norm).^2);
    %     %     t = fminbnd(polT,0,1e2)
    %     if fmin>1e-2
    %         fmin
    %         t=0:0.0001:10;
    %         polVec = zeros(length(t),1);
    %         for i=1:length(t)
    %             polVec(i) = polT(t(i));
    %         end
    %         figure;semilogy(t,real(polVec))
    %         error('Prox failed with 1e-4 tolerance')
    %
    %     end
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    % Plug-in norm of p ==t
    DeltaP = max(t-kappa,0) + 1/q_s*(max(t-kappa,0))^q_s;
    
    etaProx = eta - gamma*alpha+gamma*rho*DeltaP;
    yProx = y-gamma*p;
    
    caseInds = 2;
    
    % Case 3
elseif eta>gamma*alpha && y_norm>gamma*kappa && y_norm<=gamma*kappa+rho*(eta-gamma*alpha)
    
    etaProx = eta-gamma*alpha;
    yProx = (1-(gamma*kappa)/y_norm)*y;
    
    caseInds = 3;
    
    % Case 4
elseif eta>gamma*alpha && y_norm<=gamma*kappa
    
    etaProx = eta-gamma*alpha;
    yProx = zeros(n,1);
    
    caseInds = 4;
    
else
    warning('Case not covered')
    
    etaProx = eta;
    yProx = y;
    caseInds = 5;
    pause
end