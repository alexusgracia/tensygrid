% This scripts builds a descriptor state-space model and simulates it
% KauChr, 02.05.2025
%
% getting equilibrium point from current time instant
xEP=x0;
yEP=multilinearSimulationData.y(1,:).';
uEP=multilinearSimulationData.u(1,:).';

inputSimulation=[PrefData.';Qref*ones(1,length(out.tout));nonlinearSimulationData.Rload.';VgridDQ.';wg*ones(1,length(out.tout));out.converter2.signals.values(:,12).';zeros(1,length(out.tout))];
  
tSimIdx2=tSimIdx-1;
timeLinearization=out.tout(tSimIdx2);

            run extractingOperatingPoint;
            xOPprev=xEP;

            inputSimulationForLinearization=inputSimulation(:,tSimIdx2:tSimIdx2+1);
            inputTime=out.tout(tSimIdx2:tSimIdx2+1);
            multilinearSimulationDataForLinearization=dmsim(MTImodel,xOPprev,[],[],inputTime,inputSimulationForLinearization.',opt);
            yOPprev=multilinearSimulationDataForLinearization.y(1,:).';

   xyOP=[xEP;yEP];
            xyOPprev=[xOPprev;yOPprev];  

             xypOP=( xyOP-xyOPprev)./(out.tout(tSimIdx)-out.tout(tSimIdx2));

%% Build descriptor state-space model
A=JacobianFunctionA(timeLinearization,[xEP;yEP],[xypOP],uEP);
E=JacobianFunctionE(timeLinearization,[xEP;yEP],[xypOP],uEP);
B=JacobianFunctionB(timeLinearization,[xEP;yEP],[xypOP],uEP);
C=eye(MTImodel.n+MTImodel.p);
D=zeros(size(B));

descriptorModel=dss(A,B,C,D,E,0);

%% Construct proper descriptor model to use lsim
[msgIsProper,properDescriptorModel]=isproper(descriptorModel);

%% adjust input vector
inputDescriptorSimulation=inputSimulation-uEP;

timeVector=out.tout(tSimIdx):1e-4:out.tout(tSimIdx)+SimulationSetup.Descriptor.SimulationTime;%(length(out.tout(tSimIdx:end))-10)*1e-4;
inputDescriptorSimulation=interp1(out.tout(tSimIdx:end).',inputDescriptorSimulation(:,tSimIdx:end).',timeVector);

%% run lsim
[descriptorSimout.y,descriptorSimout.tsim,descriptorSimout.x]=lsim(properDescriptorModel,inputDescriptorSimulation,timeVector.');

%% Add equilibrium point to the output variables 
descriptorSimout.x=descriptorSimout.x.'+xEP;
descriptorSimout.y(:,1:MTImodel.n)=[descriptorSimout.y(:,1:MTImodel.n).'+xEP].';  
descriptorSimout.y(:,1+MTImodel.n:MTImodel.n+MTImodel.p)=[descriptorSimout.y(:,1+MTImodel.n:MTImodel.n+MTImodel.p).'+yEP].';  


