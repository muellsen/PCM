function betaT = proxL1o(b,lam,gam,weights)
% Proximity operator of ordered L1 norm (taken from TFOCS)

betaT = proxOrderedL1(b,lam*gam*weights);

% -- subroutines --
function x = proxOrderedL1(y,lambda)
    % Normalization
    lambda = lambda(:);
    y      = y(:);
    sgn    = sign(y);
    [y,idx] = sort(abs(y),'descend');
    
    % Simplify the problem
    k = find(y > lambda,1,'last');
    
    % Compute solution and re-normalize
    n = numel(y);
    x = zeros(n,1);
    
    if (~isempty(k))
        v1 = y(1:k);
        if numel(lambda) > 1
            v2 = lambda(1:k);
        else
            v2 = lambda*ones(k,1); % if lambda is a scalar, implicity make it lambda*ones(size(y))
        end
        v = proxAdaptiveL1Mex(v1,v2);
        x(idx(1:k)) = v;
    end
    
    % Restore signs
    x = sgn .* x;
end
end