function etaProj = projA(etaVec,gV,lb)
% Projection operator for (group) concomitant variables
% For each group the average concomitant is returned
n = length(etaVec);

% Number of groups
nG = max(gV);

etaProj = etaVec;

% Projection for the n concomitant variables
for i=1:nG
    currInds = (gV==i);
    etaProj(currInds) = mean(etaVec(currInds));
end

etaProj = max(etaProj,lb*ones(n,1));

end