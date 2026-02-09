function out = private_tensorlabChecker
    % Check if tensorlab toolbox is installed
    %
    % Description:
    %   This function checks whether the tensorlab toolbox is installed in your
    %   MATLAB environment.
    %
    % Input:
    %   None.
    %
    % Output:
    %   - installed (1 or 0): A binary indicator where '1' indicates that the
    %     tensorlab toolbox is installed, and '0' indicates that it is not
    %     installed

    try 
        % Code from tensorlabs documentation
        % https://www.tensorlab.net/doc/data.html
        hankelize(1:7,'full',false);

        T = struct;
        T.val = [1 2 3];
        T.ind = [1 2 3];
        T.incomplete = true;
        isvalidtensor(T);


        T = randn(3,5,7,9);
        M = tens2mat(T,[1 3],[4 2]);
        size(M);

        T = randn(3,5,7);
        n = 2;
        tens2mat(T,n);

        out = 1;
    catch
        out = 0;
    end
end