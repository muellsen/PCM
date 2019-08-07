function cluster = ckmeans(orivalues, n)

% first step is to sort
[values, ind] = sort(orivalues);
uniqVal =  unique(values);
cluster = ones(size(values));
clusters = {};
sortsize = length(values);
if n > length(uniqVal)
    disp('n too large, settings it to number of unique values')
    n = length(uniqVal);
end

% This is a direct translation of the javascript code in simple-statistics
% https://github.com/simple-statistics/simple-statistics/blob/master/src/ckmeans.js

if n ~= 1
    matrix =  zeros(n, sortsize);
    btmatrix = zeros(n, sortsize);
    for c=0:n-1
        first_cluster_mean = values(1);
        for sortedidx = max(c, 1):sortsize-1
            if c == 0
                sqdiff = power(values(sortedidx+1) - first_cluster_mean, 2);
                matrix(c+1,sortedidx+1) = matrix(c+1, sortedidx) + ((sortedidx / (sortedidx +1))*sqdiff);
                
                newsum = sortedidx*first_cluster_mean + values(sortedidx+1);
                first_cluster_mean = newsum /(sortedidx + 1);

            else
                sum_sq_dist = 0;
                meanXJ = 0;
                for j = sortedidx:-1:c
                    sum_sq_dist = sum_sq_dist + (sortedidx-j) / (sortedidx-j+1) * power(values(j+1) -meanXJ, 2);
                    meanXJ = (values(j+1) + (sortedidx-j)*meanXJ)/(sortedidx-j+1);

                    if j == sortedidx
                        matrix(c+1, sortedidx+1) = sum_sq_dist;
                        btmatrix(c+1, sortedidx+1) = j; 

                        if (j>0)
                            matrix(c+1, sortedidx+1) =  matrix(c, j);
                            
                        end
                    else
                        
                        if j == 0
                            if sum_sq_dist <= matrix(c+1, sortedidx+1)
                                matrix(c+1, sortedidx+1) = sum_sq_dist;
                                btmatrix(c+1, sortedidx+1) = j;
                            end
                        elseif sum_sq_dist + matrix(c, j) < matrix(c+1, sortedidx+1)
                            matrix(c+1, sortedidx+1) = sum_sq_dist + matrix(c, j);
                            btmatrix(c+1, sortedidx+1) = j;
                        end
                    end
                end
            end
        end 
    end
    tsize = (size(btmatrix) -1);
    clusterRight = tsize(2);
    
    for c = tsize(1):-1:0
        clusterLeft = btmatrix(c+1, clusterRight+1);
        clusters{c+1} = values(clusterLeft+1:clusterRight+1);
        cluster(ind(clusterLeft+1:clusterRight+1)) = c+1;
        if c > 0
            clusterRight = clusterLeft-1;
        end
    end
end
    
end