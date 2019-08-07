function [etaProx,yProx,caseInds] = proxgpBerHu(eta,y,rho,gamma)
% Perspective of BerHu in (n+1)-D

% Dimension of the problem
n = length(y);

% Norm of y
y_norm2 = sum(y.^2);
y_norm = sqrt(y_norm2);

% Check four cases for prox calculation

% Case 1
if eta > 0 && y_norm < gamma
     etaProx = eta;
     yProx = zeros(n,1);
        
     caseInds = 1;
    
% Case 2
elseif y_norm < gamma + rho*eta && y_norm > gamma % eta>0
        
        etaProx = eta;
        yProx = y-gamma*y./y_norm;
        
        caseInds = 2;
    
% Case 3
elseif 2*gamma*eta + rho*(y_norm2-gamma^2)>0 && y_norm > gamma+eta*rho
        
        % Prefactors in polynomial
        fac3 = 1;
        fac2 = 0;
        fac1 = 2/(gamma*rho^2)*(gamma+rho*eta)-1;
        fac0 = -2/(gamma*rho^2)*y_norm;
        rootfacs = [fac3 fac2 fac1 fac0];
        cubroots = roots(rootfacs);
        
        ind = intersect(find(imag(cubroots)==0),find(real(cubroots)>=0));
        
        % Show Cubic roots
        % cubroots
        
        t = cubroots(ind);        
        
        p_denom = gamma+rho*(eta+(gamma*rho/2)*(t^2-1));
        
        pp = 1/p_denom * y;
        
        if norm(pp)<1
            error('Wrong p');
        end
        
        % Return prox 
        etaProx = eta+gamma*((rho/2)*(t^2-1));
        yProx = y-gamma*pp;
        
                
        caseInds = 3;    
    
    % Case 4
elseif 2*gamma*eta + rho*(y_norm^2-gamma^2)<0 && eta<0
        
        etaProx = 0;
        yProx = zeros(n,1);
        
        caseInds = 4;
        
    else
        warning('Case not covered')
        etaProx = eta;
        yProx = y;
        caseInds = 5;
    end
end

