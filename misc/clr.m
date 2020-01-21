function clrData = clr(data,pscount)
% Centered log-ratio transform
% Input: data pxn matrix (count or compositional data)

% Add standard pseudocount if not given
if nargin<2
    pscnt = 1;
else
    pscnt = pscount;
end

% Number of variables p, Number of samples n
[p,n] = size(data);

% Column sums
sumData = sum(data);

% Check whether data is compositional
if sum(sumData(:)>1)>0
    
    % Add pseudo count 
    data = data+pscnt;

    % Normalize by the sum
    temp = data./repmat(sumData,p,1);
else
    temp = data;
end

% Calculate geometric mean
gtemp = geomean(temp);

clrData = log(temp./repmat(gtemp,p,1));


