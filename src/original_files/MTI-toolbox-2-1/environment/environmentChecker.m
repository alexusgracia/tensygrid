function systemInfo = environmentChecker(verbose, varargin)
    % Check the environment whether all dependencies are installed
    %
    % Description:
    %   This function checks whether the Tensor Toolbox is installed in your
    %   MATLAB environment.
    %
    % Input:
    %   verbose (true or false): Enable logs on the console 
    %   varargin: Optional string array with names of dependencies to check for
    %
    % Output:
    %   nargout (0): nothing
    %   nargout (1): Array of logical values indicating whether the searched addons are found in the system (1) or not found (0).
    %   The indices correspond to the indices of the optional input parameter 'array', so the first element of the return array
    %   indicates whether the first string in varargin is installed in the system or not, and so on for all subsequent positions.
    %
    % Example:
    %  a = environmentChecker(false, ["Tensorlab", "Tensor_Toolbox", "Global_Optimization_Toolbox"]);
    %  # no logs in console
    %  If Tensorlab, Tensor_Toolbox and Global_Optimization_Toolbox are installed so a is now [1, 1, 1]
    %
    %  Only Tensor_Toolbox is installed: so a is now [0, 1, 0]
    
    if nargout == 1 && nargin < 2
        fprintf(['To get a return value you have to fill the optional parameter'...
           'with a string array with the names of the searched dependencies\n'])
        fprintf('Please rerun the function and try one or more of the following input parameters: \n')
        disp(getPossibleAddonNames())
        systemInfo = -1;
        return
    end

    persistent systemInformation;


    if isempty(systemInformation) || verbose == true
        systemInformation = struct();
    end

    

    % Define persistent variable to store whether environmentChecker has been called before
    % Prevent to call this function n times, only one time via runtime is
    % required
    persistent isEnvironmentLogged
    % To reset the persistent variable enter in console "clear environmentChecker"
    
    % Check if environmentChecker has already been called
    if isempty(isEnvironmentLogged) || verbose == true
        isEnvironmentLogged = true; % Set persistent variable to true if this is the first call
    else
        if nargout == 1
            names = varargin{1};
            systemInfo = createResult(systemInformation, names);
        end
        return; % Exit function if environmentChecker has already been called
    end

    wikiPage = 'https://gitlab.cc-asp.fraunhofer.de/haw-iles/dev-ops/mti-toolbox/mti-toolbox/-/wikis/home';

    % path of json-settings-file
    filepath = 'environment/Log settings.json';

    % read json
    jsonStr = fileread(filepath);
    jsonData = jsondecode(jsonStr);

    % iterate over all addons
    for i = 1:numel(jsonData.Addon)

        % load actuall addon
        addon = jsonData.Addon(i);
        
        try
            functionToExcecute = addon.Function;
            checker = str2func(functionToExcecute);
            isPresent = checker();
        catch
            fprintf(2, 'Warning: function %s was not found from the addon %s...\n', functionToExcecute, addon.Name);
            continue
        end
       
        fieldName = strrep(addon.Name, ' ', '_');
        systemInformation.(fieldName) = isPresent;
        if verbose && addon.LogLevel == 1
            if isPresent
                disp(['Addon Found: ', addon.Name]);
            else
                fprintf(2, 'Warning: Addon %s was not found... please install it.\n', addon.Name);
				fprintf('For more information please visit the Wiki at:\n');
                fprintf('<a href="%s">%s</a>\n\n', wikiPage, wikiPage);
            end
        elseif ~verbose && addon.LogLevel == 1
            if ~isPresent
                %TODO remove duplicate code, optimize the if
                fprintf(2, 'Warning: Addon %s was not found... please install it.\n', addon.Name);
				fprintf('For more information please visit the Wiki at:\n');
                fprintf('<a href="%s">%s</a>\n\n', wikiPage, wikiPage);
            end
        end
        


    end


   %Check Matlab, after setting up CICD with licence server this file
   % and settings.json need a update
    if jsonData.MATLAB.LogLevel == 1
         if verbose 
            releaseInfo = matlabRelease;
            version = releaseInfo.Release;
            fprintf('Your Matlab version: %s matches the requirements for the mti-toolbox\n\n', version);
         end
        if isMATLABReleaseOlderThan(jsonData.MATLAB.MinimumRequiredMatlabVersion,"release") == 1 
            fprintf(2, 'Warning: Your Matlab version is outdated. Please update to version %s \n\n', jsonData.MATLAB.MinimumRequiredMarlabVersion);
        end
    end

    %Check OfficalAddon
    availableAddons = matlab.addons.installedAddons;
    for i = 1:numel(jsonData.OfficialAddon)
        
    
        %find index of selected addon in all installed addons
        addonIndex = find(strcmp(availableAddons.Name, jsonData.OfficialAddon(i).Name));
        if ~isempty(addonIndex)
            installedVersionNumber = availableAddons.Version(addonIndex);

            % check required version number
            if str2double(installedVersionNumber) >= ...
                str2double(jsonData.OfficialAddon(i).MinimumRequiredAddonVersion)
           
                if verbose && jsonData.OfficialAddon(i).LogLevel == 1
                    fprintf("Addon: %s. Version: %s\n\n",jsonData.OfficialAddon(i).Name, installedVersionNumber);
                end
                fieldName = strrep(jsonData.OfficialAddon(i).Name, ' ', '_');
                systemInformation.(fieldName) = true;
                % no verbose no print
            else
                if verbose && jsonData.OfficialAddon(i).LogLevel == 1
                    % OfficalAddon is outdated output
                    fprintf(2, ['Warning: Addon %s is installed. ' ...
                    'But your version is outdated. '...
                    'Please update to version %s \n\n'], jsonData.OfficialAddon(i).Name,...
                    jsonData.OfficialAddon(i).MinimumRequiredAddonVersion);
                    fprintf('For more information please visit the Wiki at:\n');
                    fprintf('<a href="%s">%s</a>\n\n', wikiPage, wikiPage);
                end
            end
        else
            if verbose && jsonData.OfficialAddon(i).LogLevel == 1
                % OfficalAddon is not found output
                fprintf(2, ['Warning: Addon: %s is not installed.' ...
                    ' Please install it \n'], jsonData.OfficialAddon(i).Name);
                fprintf('For more information please visit the Wiki at:\n');
                fprintf('<a href="%s">%s</a>\n\n', wikiPage, wikiPage);
            end
            fieldName = strrep(jsonData.OfficialAddon(i).Name, ' ', '_');
            systemInformation.(fieldName) = false;
        end
    end
    if nargout == 1
       names = varargin{1};
       systemInfo = createResult(systemInformation, names);
    end
end

function out = createResult(systemInformation, names)
       %Helper function that generates the result as a logical array based on the provided dependency names and system information

        if ~iscellstr(names) && ~isstring(names)
            fprintf('The optional array must be a cell array of strings.');
        end
        output = [];
        for i = 1:length(names)
             searchKey = names{i};
             if isfield(systemInformation, searchKey)
                value = systemInformation.(searchKey);
                output = [output, value]; % TODO here optimization
             else
                 fprintf('String %s, not included in the system information\n', searchKey);
                 fprintf(['Please try the following input parameters: \n'])
                 disp(getPossibleAddonNames())
             end
        end
        out = output;
end

function out = getPossibleAddonNames()
    %Helper function which returns all valid dependencies from the json file

    % path of json-settings-file
    filepath = 'environment/Log settings.json';

    % read json
    jsonStr = fileread(filepath);
    jsonData = jsondecode(jsonStr);

    % Preallocate the 'output' array for better performance
    output = cell(1, numel(jsonData.Addon) + numel(jsonData.OfficialAddon));

    % Concatenate addon names from 'jsonData.Addon'
    addonNames1 = {jsonData.Addon.Name};
    addonNames1 = strrep(addonNames1, ' ', '_');
    output(1:numel(addonNames1)) = addonNames1;

    % Concatenate addon names from 'jsonData.OfficialAddon'
    addonNames2 = {jsonData.OfficialAddon.Name};
    addonNames2 = strrep(addonNames2, ' ', '_');
    output(numel(addonNames1) + 1:end) = addonNames2;

    out = output;
end
