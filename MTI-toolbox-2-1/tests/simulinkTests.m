function tests = SimulinkTest
    tests = functiontests(localfunctions);
end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
function simulinkVersionTest(testCase)

    parentFolderPath = fileparts(fileparts(mfilename('fullpath')));
    models = dir(fullfile(parentFolderPath,filesep,"**/*.slx"));
    
    for i = 1:numel(models)
        modelFileFull = strsplit(models(i).name,'.');
        modelFile = modelFileFull{1};
    
        if ~bdIsLoaded(modelFile)
            load_system(modelFile);
        end
        
        models(i).version = get_param(modelFile,"VersionLoaded");
        
        verifyGreaterThan(testCase,models(i).version,10.5); % Oldest model currently 10.6

        close_system(modelFile);
    end


end