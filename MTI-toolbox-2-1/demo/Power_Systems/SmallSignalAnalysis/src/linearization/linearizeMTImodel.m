% Linearizing MTI model with initialization from the Simulink model
% KauChr, 25.06.2025    

% linearizing iMTI model     
            
            tic 
            run extractingOperatingPoint;
            
            xOP=xEP;
            
            inputSimulationForLinearization=inputSimulation(:,tSimIdx:end);
            inputTime=out.tout(tSimIdx:end);
            multilinearSimulationDataForLinearization=dmsim(MTImodel,xOP,[],[],inputTime(1:2),inputSimulationForLinearization(:,1:2).',opt); %%% NEEDS TO BE ADJUSTED 
            yOP=multilinearSimulationDataForLinearization.y(1,:).';

            tSimIdx2=tSimIdx-1;
            timeLinearization=out.tout(tSimIdx2);

            run extractingOperatingPoint;
            xOPprev=xEP;
            
            inputSimulationForLinearization=inputSimulation(:,tSimIdx2:end);
            inputTime=out.tout(tSimIdx2:end);
            multilinearSimulationDataForLinearization=dmsim(MTImodel,xOPprev,[],[],inputTime(1:2),inputSimulationForLinearization(:,1:2).',opt);
            yOPprev=multilinearSimulationDataForLinearization.y(1,:).';
            
            multilinearResults.GetOperatingPointTime=toc;

            
            xyOP=[xOP;yOP];
            xyOPprev=[xOPprev;yOPprev];  
            xypOP=( xyOP-xyOPprev)./(out.tout(tSimIdx)-out.tout(tSimIdx2));
            uOP=inputSimulation(:,tSimIdx);
            
            tic
            descriptorStateSpaceModel=linearize(MTImodel,xyOPprev(1:MTImodel.n),xOP,uOP,yOP,[]);
            timeLinearizationMTImodel=toc;

          
%% 
[rightEigenvectors,eigenvaluesDSS,leftEigenvectors]=eig(JacobianFunctionA(timeLinearization,[xyOP],[xypOP],uOP),JacobianFunctionE(timeLinearization,[xyOP],[xypOP],uOP));



%%
eigenvaluesDSS=eig(JacobianFunctionA(timeLinearization,[xyOP],[xypOP],uOP),JacobianFunctionE(timeLinearization,[xyOP],[xypOP],uOP));
indexNonInfinityEigenvalues=find(abs((eigenvaluesDSS))~=Inf);



