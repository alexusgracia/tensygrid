function [U,phi,rank,maxrank] = fulltens2Cpn(F,normtype,decompose, reduce)
    %FULL2CPN1 Convert a full tensor to maxrank CPN1
    %   Convert fulltensor to CPN1 with no losses. Suitable for small
    %   systems and examples. Size will grow exponentially with increasing
    %   dimensions of F.
    %  F must be a tensor of 2^(n+m) x n dimension

    arguments
        F
        normtype char = '1'
        decompose logical = false
        reduce logical = false
    end
    
    if(decompose)
        error('Decomposition not yet implemented')
    end
    
    if(normtype ~= '1')
        error("Only normtype '1' implemented as of now.")
    end
    
    dims = ndims(F);
    n = size(F,dims);
    m = dims - n - 1;
    maxrank = 2^(n+m);
    
    
    try
        % Full phi matrix (including true zeros)
        phi = reshape(F(:),[maxrank n])';
    catch
        error('Invalid tensor dimensions. Must be   2^(n+m) x n (2x2x2... for (2^(n+m) times)');
    end
    
    
    if(reduce)
        % Create a list of indices from
        all_inds = 1:2^(n+m);
    
        % Find true zeros columns
        non_empties = any(phi);
    
        % Get indeces of non zeros cols
        fulints = all_inds(non_empties);
    
        % slim out phi for true zeros
        phi = phi(:,fulints);
    end
    % determine rank
    rank = size(phi,2);
    
    % create structure-matrix by binary counting (-'0' to convert chars to
    % double)
    U = (dec2bin(0:1:2^(n+m)-1)-'0')';
    
    % reorder rows
    U(:,:) = U(end:-1:1,:);
    
    if(reduce)
        % Slim out structure matrix according to true zero columns
        U = U(:,fulints);
    end


end

