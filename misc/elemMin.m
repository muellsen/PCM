function [nnzVecUnique,outVec] = elemMin(nnzVec,inVec)

% Input: Two vectors where the first denotes the degrees of freedom and the
% second one fits
% Output: Find the vector with minimum for each degree of freedom

nnzInds = find(nnzVec(nnzVec>0));

nnzVec = nnzVec(nnzInds);
inVec = inVec(nnzInds);

nnzVecUnique = unique(nnzVec);

nUni = length(nnzVecUnique);
outVec = zeros(nUni,1);

for i=1:nUni
    
    nnzVec==nnzVecUnique(i);
    find(nnzVec==nnzVecUnique(i));
    inVec(find(nnzVec==nnzVecUnique(i)));
    outVec(i) = min(inVec(find(nnzVec==nnzVecUnique(i))));
    
end
