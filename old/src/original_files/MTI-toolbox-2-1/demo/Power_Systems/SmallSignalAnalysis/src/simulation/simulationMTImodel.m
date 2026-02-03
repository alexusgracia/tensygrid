% initialize vectors for storing the required computation time
multilinearResults.timeComputationScenarios=zeros(length(SimulationSetup.timeVectorOperatingPointsSimulation),1);
multilinearResults.errorMetrics=zeros(SimulationSetup.totalNumberOfSimulationScenarios,9);

disp(['------- Simulation of the MTI model for the case study :',SimulationSetup.Solver.MTI])

for k=1:SimulationSetup.totalNumberOfSimulationScenarios


timeLinearization=SimulationSetup.timeVectorOperatingPointsSimulation(k);

% extract intial conditions from the Simulink model
run extractingOperatingPoint.m;
x0=xEP;

% input
switch SimulationSetup.model
    case 'GFL'
        inputSimulation=[nonlinearSimulationData.Rload(tSimIdx:end).';VgridDQ(tSimIdx:end,:).';wg*ones(1,length(out.tout(tSimIdx:end)));out.converter2.signals.values(tSimIdx:end,12).';zeros(1,length(out.tout(tSimIdx:end)))];
    case 'GFM'
        inputSimulation=[nonlinearSimulationData.Rload(tSimIdx:end).';VgridDQ(tSimIdx:end,:).';wg*ones(1,length(out.tout(tSimIdx:end)));PrefData(tSimIdx:end).';Qref*ones(1,length(out.tout(tSimIdx:end)));];
    case 'GFM_GFL'
        inputSimulation=[PrefData(tSimIdx:end).';Qref*ones(1,length(out.tout(tSimIdx:end)));nonlinearSimulationData.Rload(tSimIdx:end).';VgridDQ(tSimIdx:end,:).';wg*ones(1,length(out.tout(tSimIdx:end)));out.converter2.signals.values(tSimIdx:end,12).';zeros(1,length(out.tout(tSimIdx:end)))];
end


inputTime=out.tout(tSimIdx:end);


switch SimulationSetup.Solver.MTI

     case 'dmsim'
        opt=odeset('MaxStep',SimulationSetup.Solver.MaxStep,'RelTol',SimulationSetup.Solver.RelTol);

        tic
        multilinearSimulationData=dmsim(MTImodel,x0,[],[],inputTime,inputSimulation.',opt); %stop comuptation after 10 min
        multilinearResults.timeComputationScenarios(k)=toc;
       
        
        costhetaVSMIndex=find(contains(MTImodel.stateName,'cos(thetaVSM)'));
        error_VSM=multilinearSimulationData.x(:,costhetaVSMIndex).^2+ multilinearSimulationData.x(:,costhetaVSMIndex+1).^2-1;
         
        costhetaPLLIndex=find(contains(MTImodel.stateName,'cos(thetaPLL)'));
        error_PLL=multilinearSimulationData.x(:,costhetaPLLIndex).^2+ multilinearSimulationData.x(:,costhetaPLLIndex+1).^2-1;    
        
          
        run creatingGeneralTimeseriesNames.m;

        if strcmp(SimulationSetup.model,'GFM_GFL') && ~strcmp(SimulationSetup.Descriptor.simulate,'yes')
              run plottingWithIEEEsettings.m;
        end
           
       if strcmp(SimulationSetup.Descriptor.simulate,'yes')
            run simulationDescriptorModel.m;
            run plottingDescriptorWithIEEEsettings.m; 
       end 

      
end


 disp(['Finished ',num2str(k),' simulation ouf of ',num2str(SimulationSetup.totalNumberOfSimulationScenarios),' requiring ',num2str( multilinearResults.timeComputationScenarios(k)),' s for computation.'])
end