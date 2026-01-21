%% Simulate model in abc-coordinates
% KauChr, 21.07.2025

out=sim('GFM_GFL_referenceModel_abcCoordinates');

%%
% networkDataABC=squeeze(out.networkabc.Data);
% VgridABC=networkDataABC(:,7:8);

testDataABC=out.testDataABC.signals.values;
frequencyConvVSMabc=out.omega_VSM_dq.Data/(2*pi);
angleConvVSMabc=testDataABC(:,6).';
angleGridabc=testDataABC(:,5).';

PrefData=squeeze(out.Pref.Data);%out.Pref.Data;
% extracting data
% nonlinearSimulationData.abc.Time=out.tout;
% nonlinearSimulationData.abc.Vloadabc=networkDataABC(:,5:6).';
% nonlinearSimulationData.abc.grid.IDQ=networkDataABC(:,1:2).';
% 
% nonlinearSimulationData.abc.converter1.Iabc=networkDataABC(:,3:4).';
% nonlinearSimulationData.abc.converter1.PQ=out.PQout.Data;
% 
% nonlinearSimulationData.abc.converter2.Pref=out.converter2.signals.values(:,12).';
% nonlinearSimulationData.abc.converter2.Iabc=networkDataabc(:,9:10).';
% nonlinearSimulationData.abc.converter2.PQ=out.converter2.signals.values(:,6:7).';


%%

%Simulink
    timeSimulink=out.tout.';
    % voltageDQpccSimulink=networkDataABC(:,5:6).';    
    voltagePCCAmplitudeSimulink=out.abcMAG.signals.values(:,2)*sqrt(2/3);
    gfmPQSimulink=out.PQoutABC.Data(:,:)./Sbase;
    % gfmCurrentSimulink=networkData(:,3:4).';
    gfmCurrentAmplitudeSimulink=out.abcMAG.signals.values(:,5)*sqrt(2/3);
    frequencySimulink=frequencyConvVSMabc;
    deltaAngleSimulink=(angleConvVSMabc-angleGridabc)*180/pi;

    
   gflPQSimulink=out.converter2abc.signals.values(:,6:7)./Sbase;
   % gflCurrentSimulink=networkData(:,9:10).';
   gflCurrentAmplitudeSimulink=out.abcMAG.signals.values(:,4)*sqrt(2/3);


   gflFrequencySimulink=out.converter2abc.signals.values(:,1)./(2*pi);
   gfldeltaAngleSimulink=(out.converter2abc.signals.values(:,14).')*180/pi;     


   %%

    nonlinearSimulationData.abc.Pgrid=out.abcPpu.signals.values(:,1);
    nonlinearSimulationData.abc.Pgfm=out.abcPpu.signals.values(:,2);
    nonlinearSimulationData.abc.Pgfl=out.abcPpu.signals.values(:,3);


    nonlinearSimulationData.abc.Qgrid=out.abcQpu.signals.values(:,1);
    nonlinearSimulationData.abc.Qgfm=out.abcQpu.signals.values(:,2);
    nonlinearSimulationData.abc.Qgfl=out.abcQpu.signals.values(:,3);