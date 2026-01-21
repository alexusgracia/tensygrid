function out = private_getPossibleAddonNames()
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

