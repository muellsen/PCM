function [selInds,selFreq,selMat] = stabsel(selfun,X,y,p,q,r,N,opts)
% 
%  Stability selection by Meinshausen & Buehlmann, 2010
%  Input: selfun: function handle for model selection
%         X: design matrix   
%         y: response
%         p: dimensionality of the regression vector
%         q: # of selected variables 
%         r: subsampling ratio
%         N: Number of subsamples 
%         opts: input options for selfun  
%  Output: selInds: selected indices with frequences above thresh, selFreq: selection frequencies of all
%  variables

% Number of samples 
n = size(X,1);

% Selection matrix
selMat = zeros(p,N);

% Subsample size
n_s = round(r*n);

% Thresholding numerical zeros
epsThresh = 1e-3;
    
for i=1:N
    disp([num2str(i/N*100), ' percent done...'])
    n_perm = randperm(n);
    currInds = n_perm(1:n_s);
    currX = X(currInds,:);
    currY = y(currInds);
    
    if exist('opts')
        outputMat = feval(selfun, currX, currY,opts);
    else
        outputMat = feval(selfun, currX, currY);
    end
    
    % Remove numerical zeros (if present)
    outputMat(abs(outputMat(:))<epsThresh)=0; % Thresholding the near-zero values

    % Matrix of non-zero values
    nnzMat = (outputMat~=0);
    
    if size(nnzMat,2)>1
        nnzSum = sum(nnzMat);
        % Find the first q features
        indQ = find(nnzSum<=q,1);
    else
        indQ = 1;
    end
    
    % Selected features and their regression values
    selVec = outputMat(:,indQ);
    
    % Store selected features in a matrix 
    selMat(:,i) = selVec;
    
end
    
% Empirical formula for threshold
thresh = 1/2+q^2/p;
if thresh>1
    thresh = 0.8;
end
    
selFreq = sum(selMat~=0,2)/N;
selInds = (selFreq>thresh);




