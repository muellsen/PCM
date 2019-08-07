function y = Persp_BerHu(x,eta,rho)

n = length(x);
y = zeros(n,1);

for i=1:n
   
    if eta == 0
        y(i) = abs(x(i));
        
    elseif abs(x(i))<eta*rho
        y(i) = abs(x(i));
    else
        y(i) = (x(i)^2+eta^2)/(2*eta);
    end
    
end

%y = y + 1/2*n*eta;


