% Linearizing reference model 
% KauChr, 06.06.2025

indexTimeLinearizationReferenceModel=find(out.tout>=timeLinearization,1);

%% Linearize Simulink model

tic
io=getlinio(mdl);
op=findop(mdl,timeLinearization);
nonlinearSimulationData.GetOperatingPointComputationTime=toc;
tic
referenceStateSpaceModel=linearize(mdl,op,io);
timeLinearizationSimulink=toc;



%% Compute eigenvalues
tic
referenceEigenvalues=eig(referenceStateSpaceModel);
nonlinearSimulationData.EigenvalueComputationTime=toc;