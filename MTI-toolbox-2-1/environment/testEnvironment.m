function testEnvironment(varargin)
    setup = environmentChecker(true, varargin{1});
    numberDependencies = size(varargin{1});
    data = ones(1, numberDependencies(2));
    errorMsg = 'Test failed, please install the following toolboxes on your system:_';

    if ~isequal(setup, data)
        for i = 1:length(setup)
            if setup(i) == 0
                % Note: MATLAB removes leading and trailing whitespace in a string.
                % The underscore (_) here represents a placeholder for the removed whitespace.
                errorMsg = strcat(errorMsg, varargin{1}(i), ';', '_');
            end
        end
        errorMsg = strrep(errorMsg, '_', ' ');% replace _ with whitespace
    end
    referenceString = 'Test failed, please install the following toolboxes on your system: ';
    if strcmp(errorMsg, referenceString)
        errorMsg = 'In your test, you were looking for dependencies that are not supported by the system. There might be a typo. The following dependencies are currently supported:'
        allNames = private_getPossibleAddonNames();
        for i = 1:length(allNames)
        % collect all possible addonnames
            errorMsg = [errorMsg, allNames{i}, '; '];
        end
    end

    assert(isequal(setup, data), errorMsg);
end