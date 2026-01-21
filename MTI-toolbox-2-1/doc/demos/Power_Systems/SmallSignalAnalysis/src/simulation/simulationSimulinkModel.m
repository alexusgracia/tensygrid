%%
% Simulating the reference model
% KauChr, 06.06.2025
%% choose model
switch SimulationSetup.model
    case 'GFL'
    mdl='GFL_referenceModel';
    case 'GFM'
    mdl='GFM_referenceModel';
    case 'GFM_GFL'
    %mdl='GFM_GFL_referenceModel_20250610_2023b_sorted';
    mdl='GFM_GFL_referenceModel';
end

%% no Matlab functions in Simulink model
try
tic
out=sim(mdl);
nonlinearSimulationData.SimulationComputation=toc;
catch ME
    disp('Error during simulation:');
    disp(ME.message);

    % Extract time from error message if present
    timePattern = 'at time ([\d\.eE+-]+)';
    tokens = regexp(ME.message, timePattern, 'tokens');
    if ~isempty(tokens)
        errorTime = str2double(tokens{1}{1});
        fprintf('Non-finite derivative detected at time: %g\n', errorTime);
        timeSimulationEnd=errorTime-5e-3; % Adjusting simulation end time to avoid non-finite derivative 
    end

    
    try 
    % Retry simulation with adjusted end time
    disp('Retrying simulation with adjusted end time...');
        tic
        out=sim(mdl,'StopTime',num2str(timeSimulationEnd));
        nonlinearSimulationData.SimulationComputation=toc;
    catch ME2
        disp(ME2.message);
        error('Error during adjusted simulation:');        
    end
    % return;
 end

%%
networkData=squeeze(out.network.Data);
networkData=networkData.';
VgridDQ=networkData(:,7:8);

testData=out.testData.signals.values;
frequencyConvVSM=out.omega_VSM_dq.Data/(2*pi);
angleConvVSM=testData(:,6).';
angleGrid=testData(:,5).';

PrefData=squeeze(out.Pref.Data);%out.Pref.Data;
% extracting data
nonlinearSimulationData.Time=out.tout;
nonlinearSimulationData.VloadDQ=networkData(:,5:6).';
nonlinearSimulationData.grid.IDQ=networkData(:,1:2).';
nonlinearSimulationData.Rload=out.Rload.Data;
nonlinearSimulationData.converter1.IDQ=networkData(:,3:4).';
nonlinearSimulationData.converter1.PQ=out.PQout.Data;

nonlinearSimulationData.converter2.Pref=out.converter2.signals.values(:,12).';
nonlinearSimulationData.converter2.IDQ=networkData(:,9:10).';
nonlinearSimulationData.converter2.PQ=out.converter2.signals.values(:,6:7).';