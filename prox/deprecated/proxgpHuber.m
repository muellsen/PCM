function [etaProx,yProx,caseInds] = proxgpHuber(eta,y,rho,gamma)
% Perspective of Huber in (p+1) dimensions

% Dimension of the problem
n = length(y);

% Norm of y
y_norm2 = sum(y.^2);
y_norm = sqrt(y_norm2);

% Check four cases for prox calculation

% Case 1
if (eta+y_norm^2/(2*gamma)<=0) && (y_norm <= gamma*rho)
    etaProx=0;
    yProx=zeros(n,1);
    
    caseInds = 1;
    
    % Case 2
elseif (eta <= -gamma*rho^2/2) && (y_norm>gamma*rho)
    
    etaProx = 0;
    yProx = y-gamma*rho*y./y_norm;
    caseInds = 2;
    
    % Case 3
elseif eta > -gamma*rho^2/2 && y_norm > rho*eta + gamma*rho*(1+rho^2/2);% + gamma*rho*(1+rho/2)
    
    etaProx = eta + gamma*rho^2/2;
    yProx = y  - y./y_norm*gamma*rho;
    
    caseInds = 3;
    
    % Case 4
elseif eta > -gamma*rho^2/2 && y_norm <= rho*eta + gamma*rho*(1+rho^2/2);
    
    alpha = 2;
    
    if 4*gamma*eta + 2*y_norm2 > 0
        
        % Compute mu
        mu = (4/alpha^2)*sqrt(y_norm2./(gamma^2) + (4/(27*alpha^2)*((alpha*eta)/gamma + 2)^3));
        
        p_s1 = gamma + (alpha*eta)/2;
        
        rt3term = 4*y_norm/(alpha^2*gamma);
        
        if imag(mu)==0
            p_denom = p_s1 + alpha^2*gamma/8*(sign(rt3term + mu)*abs(rt3term + mu)^(1/3) ...
                + sign(rt3term - mu)*abs(rt3term - mu)^(1/3)).^2;
        else
            p_denom = p_s1 + alpha^2*gamma/8*((rt3term + mu)^(1/3) ...
                + (rt3term - mu)^(1/3)).^2;
            warning('Imaginary part is not zero')
        end
        
        p = y/p_denom;
        p_norm = sqrt(sum(p.^2));
        
        % Test polynomial equation
        %t = y_norm/p_denom
        %res = t^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*t - 8*y_norm/(alpha^2*gamma)
        
        etaProx = eta + 1/2*gamma*(alpha*p_norm^2/2);
        yProx = y-gamma*p;
        %[etaProx,yProx] = proxgTREX(eta,y,0,0,2,gamma);
        
    else
        
        etaProx = 0;
        yProx = zeros(n,1);
        
    end
    %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
    %res2 = p^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*p - 8*y_norm/(alpha^2*gamma)
    caseInds = 4;
    
else
    warning('Case not covered')
    
    etaProx = eta;
    yProx = y;
    caseInds = 5;
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TREX prox
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [etaProx,yProx] = proxgTREX(eta,y,v,Yvec,alpha,gamma)

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);

% Shift from TREX
xTYvec = v'*Yvec;

% Shifted eta
etaPrime = eta-xTYvec;

% Check whether prox calculation is needed
if (4*gamma*etaPrime+alpha*yDiffNorm2)>0
    
    yDiffNorm = sqrt(yDiffNorm2);
    
    % Compute mu
    mu = (4/alpha^2)*sqrt(yDiffNorm2/(gamma.^2) + 32/(27*alpha^2)*((alpha*etaPrime)/(2*gamma)+1)^3);
    
    % Compute root explicitly
    p_s1 = gamma + (alpha*etaPrime)/2;
    
    rt3term = 4*yDiffNorm/(alpha^2*gamma);
    
    p_denom = p_s1 + ((alpha^2*gamma)/8)*(sign(rt3term + mu)*(abs(rt3term + mu))^(1/3) ...
        + sign(rt3term - mu)*(abs(rt3term - mu))^(1/3)).^2;
    %yDiffNorm
    p = yDiff./(p_denom);
    
    % Root checking
    %t = yDiffNorm/p_denom;
    %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
    %pause
    
    % Compute proximity
    etaProx = eta + alpha*gamma*sum((p).^2)/4;
    yProx = y-gamma*p;
    
else
    
    etaProx = xTYvec;
    yProx = Yvec;
    
end
end
