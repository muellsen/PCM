function [etaProx,yProx,caseInds] = proxg2Lq(eta,y,kappa,alpha,gamma)
% Perspective of L2 norm in (p+1) dimensions + shift
% eta: concomitant scale
% y: location 
% alpha: constant for linear operator (scale)
% kappa: constant in denominator
% gamma: scaling for prox

% Dimension of the problem
n = length(y);

% Prefactor
rho = (kappa/2);

% Norms of y
y_norm = sqrt(sum(y.^2));
y_normq = y_norm.^2;

% Default
caseInds = 2;

if y_norm==0
    etaProx = 0;
    yProx = zeros(n,1);
    return;
end

% Check cases for prox calculation
if (2*gamma*eta + rho*y_normq) > (2*gamma*alpha)
    
    % Simplified polynomial to be solved with Cardano
    % t^(3) + (2*(eta-gamma*alpha)/(gamma*rho) + q_s/(rho^2)).*t - 2*y_norm/(gamma*rho^2)==0;
    
    % Typical names for factors in Cardano's formula p amd q;
    p_c = (2*(eta-gamma*alpha)/(gamma*rho) + 2/(rho^2));
    q_c = 2*y_norm/(gamma*rho^2);
    
%     tic
%     % Simplifed polynomial
%     % t^3 + p_c*t == q_c
%     i_im = sqrt(-1);
%     
%     Q = 1/3*p_c;
%     R = 1/2*q_c;
%     D = Q^3 + R^2;
%     sD = sqrt(D);
%     
%     S = (R+sD)^(1/3);
%     T = sign(R-sD)*(abs(R-sD))^(1/3); %Modification by CLM; check before final release
%     B = S+T;
%     A = S-T;
%     
%     tVec = zeros(3,1);
%     
%     tVec(1) = B;
%     tVec(2) = -1/2*(B)+1/2*i_im*sqrt(3)*A;
%     tVec(3) = -1/2*(B)-1/2*i_im*sqrt(3)*A;
%     
%     tVec    
%     
%     t = max(tVec(abs(imag(tVec))<1e-3))
%     time1 = toc
%     
%     if isempty(t)
%         tVec;
%         t1 = max(real(tVec));
%         polR = @(t) real(t^3 + p_c.*t - q_c);
%      
% 	    t = max(0,fzero(polR,y_norm));
%         warning('Cardano root failed; using fzero to fix it')
%     end
        
    %tic
    % Companion matrix
    M = [0 0 q_c ;1 0 -p_c;0 1 0];
    qroots = eig(M);
    %time2 = toc
    
    % Only largest real root
    t = max(qroots(abs(imag(qroots))<1e-3));
    
    %t
    %pause
    
    % Simplified prox computation
    p = (y/y_norm) * t;
    
    etaProx = eta + gamma*(rho*(t^2)/2-alpha);
    
    yProx = y-gamma*p;
    
    caseInds = 1;
   
else
    etaProx = 0;
    yProx = zeros(n,1);
end



