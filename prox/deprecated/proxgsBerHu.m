function [etaProx,yProx,caseInds] = proxgsBerHu(eta,y,rho,gamma)
% Perspective of BerHu in 2D

% Dimension of the problem
p = length(y);
caseInds = zeros(p,1);

yProx = y;
etaProx = eta;

for i=1:p
    y_i = y(i);
    
    % Check four cases for prox calculation
    yabs = abs(y_i);
    
    %%
    if eta > 0 && yabs <= gamma
        
        etaProx=eta;
        yProx(i)=0;
        
        caseInds(i) = 1;
        
        %%
    elseif yabs <= gamma + rho*eta && yabs > gamma % eta>0
        
        etaProx = eta;
        yProx(i) = y_i-gamma*sign(y_i);
        
        caseInds(i) = 2;
             
        %%
    elseif 2*gamma*eta + rho*(yabs^2-gamma^2)>0 && yabs > gamma+eta*rho
        
       
        fac3 = 1;
        fac2 = 0;
        fac1 = 2/(gamma*rho^2)*(gamma+rho*eta)-1;
        fac0 = -2/(gamma*rho^2)*yabs;
        rootfacs = [fac3 fac2 fac1 fac0];
        cubroots = roots(rootfacs);
   
        ind = intersect(find(imag(cubroots)==0),find(real(cubroots)>=0))
  
        % Plot cubic root
        cubroots
        
        t = cubroots(ind)
        
        p_denom = gamma+rho*(eta+(gamma*rho/2)*(t^2-1));
        
        pp = 1/p_denom * y_i
        
        if abs(pp)<1
            error('Wrong p');
        end
        
        % Test polynomial equation
        
        etaProx = eta+gamma*((rho/2)*(pp^2-1));
        yProx(i) = y_i-gamma*pp;
        
        caseInds(i) = 3    
        
        
    elseif 2*gamma*eta + rho*max(yabs^2-gamma^2,0)<=0
        
        etaProx = 0;
        yProx(i) = 0;
        
        caseInds(i) = 4;
        
    else
        warning('Case not covered')
        etaProx = eta
        yProx(i) = y_i
        caseInds(i) = 5;
    end
        
end




