function etaProj = projC(etaVec,gV,lb)
% Projection operator for (group) concomitant variables
% with unknown grouping. The grouping is determined by 1D k-means

n = length(etaVec);

% Number of groups
nG = max(gV);
%nG1 = nG;
%gVk = gV;

% Determine optimal group
% sumVec = zeros(5,1);
% for i=1:5
%     [~,~,sumd] = kmeans(etaVec,i);
%     sumVec(i) = sum(sumd);
% end
% figure(1234)
% plot(sumVec)
% pause

% if length(unique(etaVec))>1
%     [gVk,~,~,nG1]=kmeans_opt(etaVec,2);
%     nG1
% else
%     nG1 = nG;
%     gVk = gV;
% end

%if length(unique(etaVec))>1
%    nG1 = 2;
%    gVk = ckmeans(etaVec,nG1);
    % gVk';
%else
    nG1 = nG;
    gVk = gV;
%end
    
etaProj = etaVec;

% Projection for the n concomitant variables
for i=1:nG1
    currInds = (gVk==i);
    etaProj(currInds) = mean(etaVec(currInds));
end

etaProj = max(etaProj,lb*ones(n,1));
%figure(123456)
%plot(etaProj)
%hold off
%drawnow
end