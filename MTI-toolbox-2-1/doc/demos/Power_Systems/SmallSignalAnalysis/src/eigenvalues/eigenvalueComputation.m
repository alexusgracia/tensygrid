if strcmp(SimulationSetup.model,'GFM_GFL') && strcmp(SimulationSetup.ModelABCcoordinates.simulate,'yes')
    run simulationSimulinkModel.m;
end



 switch SimulationSetup.model
                case 'GFL'
                inputSimulation=[nonlinearSimulationData.Rload.';VgridDQ.';wg*ones(1,length(out.tout));out.converter2.signals.values(:,12).';zeros(1,length(out.tout))];
                case 'GFM'
                    inputSimulation=[nonlinearSimulationData.Rload.';VgridDQ.';wg*ones(1,length(out.tout));PrefData.';Qref*ones(1,length(out.tout));];
                case 'GFM_GFL'
                inputSimulation=[PrefData.';Qref*ones(1,length(out.tout));nonlinearSimulationData.Rload.';VgridDQ.';wg*ones(1,length(out.tout));out.converter2.signals.values(:,12).';zeros(1,length(out.tout))];
                
 end


switch SimulationSetup.eigenvalueAnalysis

    case 'single'
            timeLinearization=SimulationSetup.timeVectorOperatingPointsEigenvalues(1);
            
            % reference: linearizing nonlinear Simulink model  
            run linearizeSimulinkModel.m;
            
            % linearize MTI model
            run linearizeMTImodel.m;

            % Compute eigenvalues
            eigenvaluesDSSmtiToolbox=eig(descriptorStateSpaceModel);

            tic 
            E=JacobianFunctionE(timeLinearization,xyOP,xypOP,uOP);
            A=JacobianFunctionA(timeLinearization,xyOP,xypOP,uOP);
            multilinearResults.LinearizationComputationTime=toc;

            tic
            eigenvaluesDSS=eig(A,E);
            multilinearResults.GeneralizedEigenvaluesComputationTime=toc;    
            indexNonInfinityEigenvalues=find(abs(eigenvaluesDSS)~=Inf); %only non-infinity eigenvalues
            eigenvaluesDSS=eigenvaluesDSS(indexNonInfinityEigenvalues);

            % creating a vector the same number of rows as eigenvaluesDSS
            eigenvaluesLSSlong=nan(length(eigenvaluesDSS),1);          
            eigenvaluesLSSlong(1:length(referenceEigenvalues))=referenceEigenvalues;

            % creating a comparison table and displaying it 
           
            comparisonEigenvalueTable=table(eigenvaluesLSSlong,eigenvaluesDSSmtiToolbox,'Variablenames',{['Linear SS t=',num2str(out.tout(tSimIdx)),' s'],['Descriptor SS t=',num2str(out.tout(tSimIdx)),' s']});
            disp(comparisonEigenvalueTable)

            
            
    case 'trajectory'

       trajectoryGeneralizedEigenvalues=nan(MTImodel.n,length(SimulationSetup.timeVectorEigenvalueTrajectory));
        trajectoryReferenceEigenvalues=nan(MTImodel.n,length(SimulationSetup.timeVectorEigenvalueTrajectory));


        nonlinearSimulationData.computationTimeLinearization=zeros(length(SimulationSetup.timeVectorEigenvalueTrajectory),1);
        multilinearResults.computationTimeMTIlinearization=zeros(length(SimulationSetup.timeVectorEigenvalueTrajectory),1);

       
       
         for m=1:length(SimulationSetup.timeVectorEigenvalueTrajectory)   
           
              timeLinearization=SimulationSetup.timeVectorEigenvalueTrajectory(m);

            % reference: linearizing nonlinear Simulink model  
            run linearizeSimulinkModel.m;
             nonlinearSimulationData.computationTimeLinearization(m)=timeLinearizationSimulink;

            % linearize MTI model                 
            run linearizeMTImodel.m;
            multilinearResults.computationTimeMTIlinearization(m)=timeLinearizationMTImodel;

            % Compute eigenvalues
            tic 
            E=JacobianFunctionE(timeLinearization,xyOP,xypOP,uOP);
            A=JacobianFunctionA(timeLinearization,xyOP,xypOP,uOP);
            multilinearResults.LinearizationComputationTime=toc;


            tic   
            eigenvaluesDSS=eig(A,E);
            multilinearResults.GeneralizedEigenvaluesComputationTime=toc;   
            indexNonInfinityEigenvalues=find(abs(eigenvaluesDSS)~=Inf); %only non-infinity eigenvalues
            eigenvaluesDSS=eigenvaluesDSS(indexNonInfinityEigenvalues);


           trajectoryGeneralizedEigenvalues(1:size(eigenvaluesDSS,1),m)=eigenvaluesDSS;

            % creating a vector the same number of rows as eigenvaluesDSS
            eigenvaluesLSSlong=nan(length(eigenvaluesDSS),1);          
            eigenvaluesLSSlong(1:length(referenceEigenvalues))=referenceEigenvalues;
            
           trajectoryReferenceEigenvalues(1:size(referenceEigenvalues,1),m)=referenceEigenvalues;


         end 
            
           run plottingTrajectoryEigenvaluesWithIEEEsettings.m;
    otherwise 
        disp('No computation of eigenvalues');

end