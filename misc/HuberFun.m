function [funVal,quadInds,linInds] = HuberFun(eta,res,rho,alpha,qPower)

% Number of data points
n = length(res);

q = qPower;
q_Star = 1/(1-1/qPower);

% Absolute value of residual
res_abs = abs(res);
thresh = eta.*rho.^(q_Star/q);
%find(res_abs>thresh)
linInds = find(res_abs>thresh);%,find(abs(eta)>1e-6));
quadInds = find(res_abs<thresh);%,find(abs(eta)>1e-6));
%absInds = find(abs(eta)<1e-6);


quadVals = alpha*eta(quadInds) + (res_abs(quadInds).^q)./(q.*eta(quadInds).^(q-1));
linVals = (alpha-1/q*rho^((q_Star/q))).*eta(linInds) + (rho*res_abs(linInds));
%absVals = rho*res_abs(absInds);

funVal =  sum(quadVals)+sum(linVals);%+sum(absVals);



