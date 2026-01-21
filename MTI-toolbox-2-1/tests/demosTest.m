%function tests = demosTest
%    tests = functiontests(localfunctions);
%end

% emptytest prevents warning. Mathlab need at least 1 test function
% If you implement your own tests simply remove it
%function allDemosTest(testCase)%

%    parentFolderPath = fileparts(fileparts(mfilename('fullpath')));
%    models = dir(fullfile(parentFolderPath,filesep,"**/*.mlx"));
    
%    for i = 1:numel(models)
        % This doesnnt work, as long as mlx clearr all at the beginning
        % run(models(i).name); 
%    end


%end