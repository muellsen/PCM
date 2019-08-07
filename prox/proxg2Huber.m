function [etaProx,yProx,caseInds] = proxg2Huber(eta,y,rho,alpha,gamma)
% Perspective of standard Huber in (p+1) dimensions
% eta: concomitant
% y: input
% rho: Huber's rho parameter
% alpha: constant shift
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

% Norm of y
y_norm2 = sum(y.^2);
y_norm = sqrt(y_norm2);

% Check four cases for prox calculation

% Case 1
if (y_norm<=gamma*rho) && (y_norm2 <= 2*(gamma^2)*(alpha-eta/gamma))
    
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif (y_norm>gamma*rho) && (eta <= gamma*(alpha-(rho^2)/2))
    
    etaProx = 0;
    yProx = (1-(gamma*rho)/y_norm)*y;
    caseInds = 2;
    
    % Case 3
elseif (eta > gamma*(alpha-(rho^2)/2)) && y_norm >= gamma*rho*(eta/gamma+1+(rho^2)/2-alpha)
    
    etaProx = eta + gamma*((rho^2)/2-alpha);
    yProx = (1-(gamma*rho)/y_norm)*y;
    
    caseInds = 3;
    
    % Case 4
elseif (y_norm2 > 2*(gamma^2)*(alpha-(eta/gamma))) && (y_norm < gamma*rho *(eta/gamma + 1 + (rho^2)/2 - alpha));
    
    % Standard companion matrix
    M = [[zeros(1,2);eye(2,2)],zeros(3,1)];
    
    % Explicit root finding via determinant method
    M(1,3) = 2*y_norm/gamma;
    M(2,3) = -(2.*(eta-gamma*alpha)+2*gamma)/gamma;
    qroots = eig(M);
    
    % Only largest real root
    t = max(qroots(abs(imag(qroots))<1e-3));
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    % Redundant (just for legacy)
    %if 2*gamma*eta + y_norm2 > 2*(gamma^2)*alpha 
    
    etaProx = (eta + gamma*((t.^2)/2-alpha));
    yProx = y-gamma*p;
    
    %     else
    %
    %         etaProx = 0;
    %         yProx = zeros(n,1);
    %
    %     end
    
    caseInds = 4;
    
else
    warning('Case not covered')
    
    etaProx = eta
    yProx = y
    caseInds = 5;
    pause
end
