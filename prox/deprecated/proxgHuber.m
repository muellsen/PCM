function [etaProx,yProx,caseInds] = proxgHuber(eta,y,rho,gamma)
% Perspective of Huber in 2D

% Dimension of the problem
p = length(y);
caseInds = zeros(p,1);

yProx = y;
etaProx = eta;

for i=1:p
    y_i = y(i);
    
    % Check five cases for prox calculation
    yabs = abs(y_i);
    
    % If ? + |y|^2/(2?) < 0 and |y| <= ??, then Theorem 3.1(i) yields (?, q) = (0, 0).
    
    if (eta+yabs^2/(2*gamma)<=0) && (yabs <= gamma*rho)
        etaProx=0;
        yProx(i)=0;
        
        caseInds(i) = 1;
        
        % If ? <= ???^2/2 and |y| > ??
    elseif (eta <= -gamma*rho^2/2) && (yabs>gamma*rho)
        
        etaProx = 0;
        yProx(i) = y_i-gamma*rho*sign(y_i);
        
        caseInds(i) = 2;
        
        % If ? > ???^2/2 and y > ? + ??(1 + ?/2)
    elseif eta > -gamma*rho^2/2 && y_i > rho*eta + gamma*rho*(1+rho^2/2);% + gamma*rho*(1+rho/2)
        
        etaProx = eta + gamma*rho^2/2;
        yProx(i) = y_i - gamma*rho;
        
        caseInds(i) = 3;
        % If ? > ???^2/2 and y < - ? - ??(1 + ?/2)
    elseif eta > -gamma*rho^2/2 && y_i < -rho*eta - gamma*rho*(1+rho^2/2);% - gamma*rho*(1+rho/2)
        
        etaProx = eta + gamma*rho^2/2;
        yProx(i) = y_i + gamma*rho;
        
        caseInds(i) = 4;
        %  ?/? > ??2/2 and |y| <=  ? + ??(1 + ?/2)
        
    elseif eta > -gamma*rho^2/2 && yabs <= rho*eta + gamma*rho*(1+rho^2/2);
                
        yNorm2 = sum(y_i.^2);
        yNorm = sqrt(yNorm2);
        alpha = 2;
        
        if 4*gamma*eta + 2*yNorm2 > 0
            
            % Compute mu
            mu = (4/alpha^2)*sqrt(yNorm2./(gamma^2) + (4/(27*alpha^2)*((alpha*eta)/gamma + 2)^3));
            
            p_s1 = gamma + (alpha*eta)/2;
            
            rt3term = 4*yNorm/(alpha^2*gamma);
            
            p_denom = p_s1 + alpha^2*gamma/8*(sign(rt3term + mu)*abs(rt3term + mu)^(1/3) ...
                + sign(rt3term - mu)*abs(rt3term - mu)^(1/3)).^2;
            p = y_i/p_denom;
            p_norm = sqrt(sum(p.^2));
            
            % Test polynomial equation
            %t = yNorm/p_denom
            %res = t^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*t - 8*yNorm/(alpha^2*gamma)
            
            etaProx = eta + 1/2*gamma*(alpha*p_norm^2/2);
            yProx(i) = y_i-gamma*p;
            %[etaProx,yProx] = proxgTREX(eta,y,0,0,2,gamma);
            
        else
            
            etaProx = 0;
            yProx(i) = 0;
            
        end
        %t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
        %res2 = p^3+(4*(alpha*eta+2*gamma)./(alpha^2*gamma))*p - 8*yNorm/(alpha^2*gamma)
        
        caseInds(i) = 5;
    else
        warning('Case not covered')
        
        etaProx = eta;
        yProx = y_i;
        caseInds(i) = 6;
        
    end
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
