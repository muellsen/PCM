function [etaProx,yProx,caseInds] = proxgVapnik(eta,y,eps,gamma)
% Perspective of Vapnik in (n+1) dimensions

% Dimension of the problem
n = length(y);

% Check seven cases for prox calculation
y_norm = sqrt(sum(y.^2));

if (eta+eps*y_norm<=0) && (y_norm <= gamma)
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;

elseif (eta<=-gamma*eps) && (y_norm > gamma)
    
    etaProx = 0;
    yProx = y-gamma*y./y_norm;
    
    caseInds = 2;
    
elseif eta > -gamma*eps && y_norm > eps*eta + gamma*(1+eps^2);
    
    etaProx = eta + gamma*eps;
    yProx = y - gamma*y./y_norm;
    
    caseInds = 3;
    
elseif y_norm>-eta/eps && eps*eta<= y_norm  && y_norm <=eps*eta + gamma*(1+eps^2)
    
    etaProx = (eta+eps*y_norm)/(1+eps^2);
    yProx = eps*(eta+eps*y_norm)*y/(y_norm*(1+eps^2));
    
    caseInds = 4;
    
elseif eta >= 0 && y_norm <= eps*eta;
    
    etaProx = eta;
    yProx = y;
    
    caseInds = 5;
    
else
    
    warning('Case not covered')
    
    etaProx = eta;
    yProx = y;
    caseInds = 6;
    
end



