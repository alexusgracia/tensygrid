%
% This scripts creates a new variables from the obtained simulation data to
% to allow for a generic plotting script. 
% KauChr, 01.05.2025


                     %Simulink
                    timeSimulink=out.tout.';
                    voltageDQpccSimulink=networkData(:,5:6).';    
                    voltagePCCAmplitudeSimulink=sqrt(voltageDQpccSimulink(1,:).^2+voltageDQpccSimulink(2,:).^2)./Vpeak;
                    gfmPQSimulink=out.PQout.Data(:,:)./Sbase;
                    gfmCurrentSimulink=networkData(:,3:4).';
                    gfmCurrentAmplitudeSimulink=sqrt(gfmCurrentSimulink(1,:).^2+gfmCurrentSimulink(2,:).^2)./Ibase;
                    frequencySimulink=frequencyConvVSM;
                    deltaAngleSimulink=(angleConvVSM-angleGrid)*180/pi;

                    

                   gflPQSimulink=out.converter2.signals.values(:,6:7)./Sbase;
                   gflCurrentSimulink=networkData(:,9:10).';
                   gflCurrentAmplitudeSimulink=sqrt(gflCurrentSimulink(1,:).^2+gflCurrentSimulink(2,:).^2)./Ibase;


                   gflFrequencySimulink=out.converter2.signals.values(:,1)./(2*pi);
                   gfldeltaAngleSimulink=(out.converter2.signals.values(:,13).')*180/pi;     


           switch SimulationSetup.model

               case 'GFM'
                   indexGRIDiD=find(MTImodel.stateName=='grid_ID');
                    indexGRIDangle=find(MTImodel.stateName=='grid_thetaGRID');
                    indexVDload=find(MTImodel.algebraicName=='load1_VD');

                    indexGFMiD=find(MTImodel.stateName=='converter1_ID');
                    indexGFMp=find(MTImodel.algebraicName=='converter1_Pmeas');
                    indexGFMomega=find(MTImodel.algebraicName=='converter1_omegaVSM');
                    indexGFMangle=find(MTImodel.stateName=='converter1_thetaVSM');

                    timeDMSIM=multilinearSimulationData.tsim;
                    voltageAmplitudeDMSIM=sqrt((multilinearSimulationData.y(:,indexVDload).^2).'+(multilinearSimulationData.y(:,indexVDload+1).^2).')./Vpeak; % perhaps better with sum or accumarray
                   
                    gfmCurrentAmplitudeDMSIM=sqrt((multilinearSimulationData.x(:,indexGFMiD).^2).'+(multilinearSimulationData.x(:,indexGFMiD+1).^2).')./Ibase;
                    gfmPQDMSIM=multilinearSimulationData.y(:,indexGFMp:indexGFMp+1)./Sbase;
                    gfmFrequencyDMSIM=multilinearSimulationData.y(:,indexGFMomega)/(2*pi);
                    deltaAngleDMSIM=(multilinearSimulationData.x(:,indexGFMangle)-multilinearSimulationData.x(:,indexGRIDangle))*180/pi;
               case 'GFM_GFL'
                    
                    indexGRIDiD=find(MTImodel.stateName=='grid_ID');
                    indexGRIDangle=find(MTImodel.stateName=='grid_thetaGRID');
                    indexVDload=find(MTImodel.algebraicName=='load1_VD');

                    indexGFMiD=find(MTImodel.stateName=='converter1_ID');
                    indexGFMp=find(MTImodel.algebraicName=='converter1_Pmeas');
                    indexGFMomega=find(MTImodel.algebraicName=='converter1_omegaVSM');
                    indexGFMangle=find(MTImodel.stateName=='converter1_thetaVSM');

                    indexGFLiD=find(MTImodel.stateName=='converter2_ID');
                    indexGFLp=find(MTImodel.algebraicName=='converter2_Pmeas');
                    indexGFLomega=find(MTImodel.algebraicName=='converter2_omegaPLL');
                    indexGFLangle=find(MTImodel.stateName=='converter2_thetaPLL');


                    %dmsim
                    timeDMSIM=multilinearSimulationData.tsim;
                    voltageAmplitudeDMSIM=sqrt((multilinearSimulationData.y(:,indexVDload).^2).'+(multilinearSimulationData.y(:,indexVDload+1).^2).')./Vpeak; % perhaps better with sum or accumarray
                   
                    gfmCurrentAmplitudeDMSIM=sqrt((multilinearSimulationData.x(:,indexGFMiD).^2).'+(multilinearSimulationData.x(:,indexGFMiD+1).^2).')./Ibase;
                    gfmPQDMSIM=multilinearSimulationData.y(:,indexGFMp:indexGFMp+1)./Sbase;
                    gfmFrequencyDMSIM=multilinearSimulationData.y(:,indexGFMomega)/(2*pi);
                    deltaAngleDMSIM=(multilinearSimulationData.x(:,indexGFMangle)-multilinearSimulationData.x(:,indexGRIDangle))*180/pi;

                    gflCurrentAmplitudeDMSIM=sqrt((multilinearSimulationData.x(:,indexGFLiD).^2).'+(multilinearSimulationData.x(:,indexGFLiD+1).^2).')./Ibase;
                    gflPQDMSIM=multilinearSimulationData.y(:,indexGFLp:indexGFLp+1)./Sbase;
                    gflFrequencyDMSIM=multilinearSimulationData.y(:,indexGFLomega)/(2*pi);
                    gfldeltaAngleDMSIM=(multilinearSimulationData.x(:,indexGFLangle)-multilinearSimulationData.x(:,indexGRIDangle))*180/pi;

               case 'GFL'
                    indexGRIDiD=find(MTImodel.stateName=='grid_ID');
                    indexGRIDangle=find(MTImodel.stateName=='grid_thetaGRID');
                    indexVDload=find(MTImodel.algebraicName=='load1_VD');

                    indexGFLiD=find(MTImodel.stateName=='converter2_ID');
                    indexGFLp=find(MTImodel.algebraicName=='converter2_Pmeas');
                    indexGFLomega=find(MTImodel.algebraicName=='converter2_omegaPLL');
                    indexGFLangle=find(MTImodel.algebraicName=='converter2_thetaPLL');

                   timeDMSIM=multilinearSimulationData.tsim;
                    voltageAmplitudeDMSIM=sqrt((multilinearSimulationData.y(:,indexVDload).^2).'+(multilinearSimulationData.y(:,indexVDload+1).^2).')./Vpeak; % perhaps better with sum or accumarray
                    
                     gflCurrentAmplitudeDMSIM=sqrt((multilinearSimulationData.x(:,indexGFLiD).^2).'+(multilinearSimulationData.x(:,indexGFLiD+1).^2).')./Ibase;
                    gflPQDMSIM=multilinearSimulationData.y(:,indexGFLp:indexGFLp+1)./Sbase;
                    gflFrequencyDMSIM=multilinearSimulationData.y(:,indexGFLomega)/(2*pi);
                    gfldeltaAngleDMSIM=(multilinearSimulationData.x(:,indexGFLangle)-multilinearSimulationData.x(:,indexGRIDangle))*180/pi;

           end
                    