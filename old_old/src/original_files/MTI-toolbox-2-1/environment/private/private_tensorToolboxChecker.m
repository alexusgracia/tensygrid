function out = private_tensorToolboxChecker
    % Check if tensor toolbox is installed
    %
    % Description:
    %   This function checks whether the tensor toolbox is installed in your
    %   MATLAB environment.
    %
    % Input:
    %   None.
    %
    % Output:
    %   - installed (1 or 0): A binary indicator where '1' indicates that the
    %     Tensor Toolbox is installed, and '0' indicates that it is not
    %     installed

    try
        % Code from tensor toolbox documentation
        % https://www.tensortoolbox.org/tensor_doc.html
        M = ones(4,3,2); %<-- A 4 x 3 x 2 array.
        tensor(M); %<-- Convert to a tensor object.

        tensor(M,[2 3 4]); %<-- M has 24 elements.

        tensor(rand(5,1)); %<-- Creates a 2-way tensor.

        tensor(rand(4,3,1),[4 3 1]); %<-- Creates a 3-way tensor.

        tensor; %<-- Creates an empty tensor.

        X = tenrand([3 4 2 1]); %<-- Create a 3 x 4 x 2 x 1 random tensor.
        X(1,1,1,1); %<-- Extract a single element.

        out = 1;
    catch
        out = 0;
    end

end