function x = prox_l1sum( x0, lambda, b )

    brk_pts = sort( [x0-lambda;x0+lambda], 'descend' );

    shrink  = @(x) sign(x).*max( abs(x) - lambda, 0 );
    xnu     = @(nu) shrink( x0 - nu );
    h       = @(x) sum(x) - b; % want to solve h(nu) = 0

    % Bisection
    lwrBnd       = 0;
    uprBnd       = length(brk_pts) + 1;
    iMax         = ceil( log2(length(brk_pts)) ) + 1;
    PRINT = false; % set to "true" for debugging purposes
    if PRINT
        dispp = @disp;
        printf = @fprintf;
    else
        dispp = @(varargin) 1;
        printf = @(varargin) 1;
    end
    dispp(' ');
    for i = 1:iMax
        if uprBnd - lwrBnd <= 1
            dispp('Bounds are too close; breaking');
            break;
        end
        j = round(mean([lwrBnd,uprBnd]));
        printf('j is %d (bounds were [%d,%d])\n', j, lwrBnd,uprBnd );
        if j==lwrBnd
            dispp('j==lwrBnd, so increasing');
            j = j+1;
        elseif j==uprBnd
            dispp('j==uprBnd, so increasing');
            j = j-1;
        end
        
        a   = brk_pts(j);
        x   = xnu(a);  % the prox
        p   = h(x);
        
        if p > 0
            uprBnd = j;
        elseif p < 0
            lwrBnd = j;
        end
        if PRINT
            % Don't rely on redefinition of printf,
            % since then we would still calculate find(~x)
            % which is slow
            printf('i=%2d, a = %6.3f, p = %8.3f, zeros ', i, a, p );
            if i < 100, printf('%d ', find(~x) ); end
            printf('\n');
        end
    end
    
    % Now, determine linear part, which we infer from two points.
    % If lwr/upr bounds are infinite, we take special care
    % e.g., we make a new "a" slightly lower/bigger, and use this
    % to extract linear part.
    %uprBnd

    if lwrBnd == 0
        a2 = brk_pts( uprBnd );
        a1 = a2 - 10; % arbitrary
        aBounds = [-Inf,a2];
    elseif uprBnd == length(brk_pts) + 1
        a1 = brk_pts( lwrBnd );
        a2 = a1 + 10; % arbitrary
        aBounds = [a1,Inf];
    else
        % In general case, we can infer linear part from the two break points
        a1 = brk_pts( lwrBnd );
        a2 = brk_pts( uprBnd );
        aBounds = [a1,a2];
    end
    
    % Now we have the support, find exact value
    x       = xnu( mean(aBounds) );  % to find the support
    supp    = find(x);

    sgn     = sign(x);
    nu      = ( sum(x0(supp) - lambda*sgn(supp) ) - b )/length(supp);
    
    x   = xnu( nu );

end