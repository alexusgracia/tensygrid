function out = private_sparseToolkitChecker
% Check if sparse toolkit is installed and its version
%
% Description:
%   This function checks whether the sparse toolkit is installed in your
%   MATLAB environment.
%
% Input:
%   None.
%
% Output:
%   - installed (1 or 0): A binary indicator where '1' indicates that the
%     sparse toolkit is installed, and '0' indicates that it is not
%     installed or the version is not right

    try 
        % Code from tensorlabs documentation
        % https://sites.google.com/view/sparse-grids-kit sparse grids
        % MATLAB kit
        % Test functionality
        level = 3;
        d = 4;
        knots=@(n) knots_CC(n,-1,1,'nonprob'); % Clenshaw-Curtis knots
        w = level;
        S = create_sparse_grid(d,w,knots,@lev2knots_doubling); % this function is only valid from realse 23-5 ("Robert") or above.
        reduce_sparse_grid(S);
    

        out = 1;
    catch
        out = 0;
    end
end