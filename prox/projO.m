function etaProj = projO(etaVec,gV,lb)
% Projection operator for organic Lasso
% For each group the vector on the unit simplex is returned
n = length(etaVec);

% Number of groups
nG = max(gV);

etaProj = max(etaVec,zeros(n,1));

% Projection for the n concomitant variables
for i=1:nG
    currInds = (gV==i);
    etaProj(currInds) = etaProj(currInds)./sum(etaProj(currInds));
end

end