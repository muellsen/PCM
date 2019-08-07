function [etaProx,yProx,caseInds] = proxg2BerHu(eta,y,rho,alpha,kappa,gamma)
% Perspective of standard BerHu in (p+1) dimensions
% eta: concomitant scale
% y: location 
% rho: BerHu's rho parameter
% alpha: constant in linear part
% kappa: Slope of the inner BerHu part
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

% Delta(y/gamma)

ygamnorm = sqrt(sum((y/gamma).^2));
y_norm = sqrt(sum(y.^2));

Delta_ygam = max(ygamnorm-kappa,0) + 1/2*(max(ygamnorm-kappa,0))^2;

% Check four cases for prox calculation

% Case 1
if Delta_ygam<=(alpha-eta/gamma)/rho
    
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif y_norm>=gamma*kappa+rho*(eta-gamma*alpha) && Delta_ygam > (alpha-eta/gamma)/rho
    
    % Polynomial t^3 + p_c*t - q_c = 0
    p_c = 2*(gamma+rho*(eta-gamma*alpha))/(gamma*rho^2)-1;
    q_c = 2*y_norm/(gamma*rho^2);
    
    M = [[zeros(1,2);eye(2,2)],zeros(3,1)];
    
    % Explicit root finding via determinant method
    M(1,3) = q_c;
    M(2,3) = -p_c;
    qroots = eig(M);
    
    % Only largest real root
    t = max(qroots(abs(imag(qroots))<1e-3));
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    % Plug-in norm of p==t
    DeltaP = max(t-kappa,0) + 1/2*(max(t-kappa,0))^2;
    
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