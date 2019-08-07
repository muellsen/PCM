function [etaProx,yProx,caseInds] = proxgqHuber(eta,y,rho,alpha,qPower,gamma)
% Perspective of generalized Huber in (p+1) dimensions
% eta: concomitant scale
% y: location 
% rho: Huber's rho parameter
% alpha: constant
% qPower: qth power of the inner Huber part
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

q_s = qPower/(qPower-1);

% Norm of y
y_norm = sqrt(sum(y.^2));
y_normq = y_norm^q_s;

% Check four cases for prox calculation

% Case 1
if (y_norm<=gamma*rho) && (y_normq <= q_s*(gamma^q_s)*(alpha-eta/gamma))
    
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif (y_norm>gamma*rho) && (eta <= gamma*(alpha-(rho^q_s)/q_s))
    
    etaProx = 0;
    yProx = (1-gamma*rho/y_norm)*y;
    caseInds = 2;
    
    % Case 3
elseif (eta > gamma*(alpha-(rho^q_s)/q_s)) && y_norm >= gamma*rho^(q_s-1)*(eta/gamma+rho^(2-q_s)+(rho^q_s)/q_s-alpha)
    
    etaProx = eta + gamma*((rho^q_s)/q_s-alpha);
    yProx = (1-(gamma*rho)/y_norm)*y;
    
    caseInds = 3;
    
    % Case 4
elseif (y_normq > q_s*(gamma^q_s)*(alpha-(eta/gamma))) && (y_norm < gamma*rho^(q_s-1) *(eta/gamma + rho^(2-q_s) + (rho^q_s)/q_s - alpha));
    
    
    % Real part of the polynomial (for numerical reasons only)
    polR = @(t) real(t^(2*q_s-1) + q_s/gamma*(eta-gamma*alpha)*t^(q_s-1)+q_s*t - q_s*y_norm/gamma);
    
    t = fzero(polR,y_norm);
    
    fmin = polR(t);
    
    if fmin>1e-4
        error('Prox failed with 1e-4 tolerance')
    end
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    if q_s*gamma^(q_s-1)*eta + y_normq > q_s*(gamma^q_s)*alpha
        
        etaProx = (eta + gamma*((t.^q_s)/q_s-alpha));
        yProx = y-gamma*p;
    
    else
        
        etaProx = 0;
        yProx = zeros(n,1);
    
    end
    
    caseInds = 4;
    
else
    warning('Case not covered')
    
    etaProx = eta
    yProx = y
    caseInds = 5;
    pause
end
