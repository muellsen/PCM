function [etaProx,yProx,caseInds] = proxgqLq(eta,y,kappa,alpha,qPower,gamma)
% Perspective of generalized Lq norm in (p+1) dimensions

% eta: concomitant scale
% y: location 
% alpha: constant for linear operator (scale)
% kappa: constant in denominator
% qPower: qth power of the norm
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

% Dual exponent
q_s = qPower/(qPower-1);

% Prefactor
rho = (kappa/qPower)^(q_s-1);

% Norms of y
y_norm = sqrt(sum(y.^2));
y_normq = y_norm.^q_s;

% Default
caseInds = 2;

if y_norm==0
    etaProx = 0;
    yProx = zeros(n,1);
    return;
end

% Check cases for prox calculation
if (q_s*gamma^(q_s-1)*eta + rho*y_normq) > (q_s*gamma^q_s*alpha)
    
    % Real part of the polynomial (for numerical reasons only)
    polR = @(t) real( t^(2*q_s-1) + q_s*(eta-gamma*alpha)/(gamma*rho).*t^(q_s-1) + q_s/(rho^2).*t - q_s*y_norm/(gamma*rho^2));
    
    t = fzero(polR,y_norm);
    
    fmin = polR(t);
    
    if fmin>1e-4
        error(['Prox failed with 1e-4 tolerance: min_f = ',num2str(fmin)])
    end
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    etaProx = eta + gamma*(rho*(t^q_s)/q_s-alpha);
    yProx = y-gamma*p;
    
    caseInds = 1;
else
    etaProx = 0;
    yProx = zeros(n,1);
end
    
    

