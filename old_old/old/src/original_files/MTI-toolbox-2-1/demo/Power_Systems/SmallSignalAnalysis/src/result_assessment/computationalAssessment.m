%% Assess the computational aspects of the simullation, linearization of the iMTI and the NTI model, as well as the commputation of the generalized eigenvalues and the eigenvalues 
% KauChr, 23.07.2025


multilinearResults.SimulationComputationTime=multilinearResults.timeComputationScenarios(find(multilinearResults.timeComputationScenarios~=0));

rowNames={'Simulation'};
nonlinearSimulationData.Summary=[nonlinearSimulationData.SimulationComputation];
multilinearResults.Summary=[multilinearResults.SimulationComputationTime];



if ~strcmp(SimulationSetup.eigenvalueAnalysis,'no')
    nonlinearSimulationData.totalComputationTimeLinearizationSimulink=sum(nonlinearSimulationData.computationTimeLinearization);
    nonlinearSimulationData.averageComputationTimeLinearizationSimulink=mean(nonlinearSimulationData.computationTimeLinearization);
    multilinearResults.totalComputationTimeLinearizationMTI=sum(multilinearResults.computationTimeMTIlinearization);
    multilinearResults.averageComputationTimeLinearizationMTI=mean(multilinearResults.computationTimeMTIlinearization);

    rowNames={'Simulation','Jacobian','linearize','Operating Point','Eigenvalues'};
    
    nonlinearSimulationData.Summary=[nonlinearSimulationData.SimulationComputation;nan;nonlinearSimulationData.averageComputationTimeLinearizationSimulink;nonlinearSimulationData.GetOperatingPointComputationTime;nonlinearSimulationData.EigenvalueComputationTime]*1e3;
    % Handle missing or empty JacobianTime field
    if isfield(multilinearResults, 'JacobianTime') && ~isempty(multilinearResults.JacobianTime)
        jacobianTimeValue = multilinearResults.JacobianTime;
    else
        jacobianTimeValue = nan;
    end
    multilinearResults.Summary=[multilinearResults.SimulationComputationTime;jacobianTimeValue;multilinearResults.LinearizationComputationTime;multilinearResults.GetOperatingPointTime ;multilinearResults.GeneralizedEigenvaluesComputationTime]*1e3;
end 

computationalAssessmentTable=table(nonlinearSimulationData.Summary,multilinearResults.Summary,'RowNames',rowNames,'VariableNames',{'NTI (ms)','iMTI (ms)'});
disp(computationalAssessmentTable)