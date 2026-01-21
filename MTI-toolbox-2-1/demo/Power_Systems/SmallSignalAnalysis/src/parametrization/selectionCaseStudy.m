

switch SimulationSetup.caseStudy

    case 'SCR5_XR10_GFL_inverter_Pref_steps'

           SimulationSetup.Solver.timeSimulationEnd=4; 

            
           % load 
           SimulationSetup.Load.Base=Rload;
           SimulationSetup.Load.Step.Time=10;
           SimulationSetup.Load.Step.Value=Rload;
            



           SimulationSetup.grid.voltageAmplitudeStep.Time=3.4;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92];

           SimulationSetup.grid.voltageAmplitudeStep.Time=3.4;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];%[0;(0.01*0.0174533)*1.01];;

           SimulationSetup.converter1.Pref.Time=3;
           SimulationSetup.converter1.Pref.Values=[0.001;0.001];

           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.1;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter2.Pref.Time=[0,0.1+[1  ,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,  2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9]];
           SimulationSetup.converter2.Pref.Values=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1  ,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];

           SimulationSetup.timeVectorOperatingPointsSimulation=[1.099];%[2.199]%[0.995];%[2.95;2.95;2.9;2.8];%[2.999;2.95;2.9;2.8];
%1.199


            SimulationSetup.timeVectorOperatingPointsEigenvalues=1.099;%1.205%0.995;
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.1:1.5].';
%3.4+[0:0.1:0.4].'%


    case 'SCR5_XR10_GFL_rectifier_Pref_steps' 

           SimulationSetup.Solver.timeSimulationEnd=4;  

           SimulationSetup.grid.voltageAmplitudeStep.Time=3.4;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92];

           timePrefStep=3;
           PrefStep=[0.001;0.001]%[0.01375; 0.7];

           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.1;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter2.Pref.Time=[0,0.1+[1  ,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,  2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9]];
           SimulationSetup.converter2.Pref.Values=-[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1  ,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];

    case 'SCR1_XR10'

       SimulationSetup.Solver.timeSimulationEnd=2.8; 

           % load 
           SimulationSetup.Load.Base=Rload*0.5;
           SimulationSetup.Load.Step.Time=2.5;
           SimulationSetup.Load.Step.Value=Rload;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=2.6;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.2;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=1.2;
           SimulationSetup.converter1.Pref.Values=[0.001;0.5];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,1.4];
           SimulationSetup.converter2.Pref.Values=[0,0.5];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[2.5-1e-3];

            
            SimulationSetup.timeVectorOperatingPointsEigenvalues=[2.5-1e-3];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.05:0.3].'; 

        SCR=1;                   % short circuit ratio
        XR_ratio=10;             % X/R ratio  


        Scc=SCR*Sbase;           % short circuit capacity
        Xcct=Vbase^2/Scc;         % reactance 
        Rg_t=sqrt(Xcct^2/(XR_ratio^2+1));      
        Xg_t=Rg_t*XR_ratio;   

        % Thevenin equivalent grid 
        Zg_n2=complex(Rg_t,Xg_t);
        Rg=real(Zg_n2);         % Thevenin resistance
        Lg=imag(Zg_n2)/wg;
        Zgrid_pu=sqrt(Rg^2+Lg^2)/Zbase;
        

    case 'GFLonly_GEtest'
            
           
           SimulationSetup.Solver.timeSimulationEnd=2; 

           % load 
           SimulationSetup.Load.Base=Rload;
           SimulationSetup.Load.Step.Time=3;
           SimulationSetup.Load.Step.Value=Rload*1.2;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=1.5;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 


           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];%[0;(0.01*0.0174533)*1.01];;

           SimulationSetup.converter1.Pref.Time=3;
           SimulationSetup.converter1.Pref.Values=[0.001;0.001];

           SimulationSetup.converter2.Pref.Time=[0,0.1+[1  ,1.1,1.2,1.3,1.4 ]];
           SimulationSetup.converter2.Pref.Values=[0,0.1,0.2,0.3,0.4,0.5];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[1.499];

            SimulationSetup.timeVectorOperatingPointsEigenvalues=[1.499];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.005:0.05].';

    case 'GFMonly'

           SimulationSetup.Solver.timeSimulationEnd=2; 

           % load 
           SimulationSetup.Load.Base=Rload;
           SimulationSetup.Load.Step.Time=3;
           SimulationSetup.Load.Step.Value=Rload*1.2;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=1.5;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.1;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=2;
           SimulationSetup.converter1.Pref.Values=[0.001;0.6];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,0.1+[1  ,1.1,1.2,1.3,1.4 ]];
           SimulationSetup.converter2.Pref.Values=[0,0.1,0.2,0.3,0.4,0.5];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[1.499];

            SimulationSetup.timeVectorOperatingPointsEigenvalues=[1.499];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.005:0.05].';


    case 'VoltageStep_SCR5_XR10_equalLoadDistribution_Rstep055'
          SimulationSetup.Solver.timeSimulationEnd=3; 

           % load 
           SimulationSetup.Load.Base=Rload*0.5;
           SimulationSetup.Load.Step.Time=2.5;
           SimulationSetup.Load.Step.Value=Rload*0.55;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=2;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.2;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=1.2;
           SimulationSetup.converter1.Pref.Values=[0.001;0.5];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,1.4];
           SimulationSetup.converter2.Pref.Values=[0,0.5];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[2.5-1e-3];

            
            SimulationSetup.timeVectorOperatingPointsEigenvalues=[2.499];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.005:0.05].';
            
    case 'VoltageStep_SCR5_XR10_equalLoadDistribution_Rstep'
          SimulationSetup.Solver.timeSimulationEnd=2.8; 

           % load 
           SimulationSetup.Load.Base=Rload*0.5;
           SimulationSetup.Load.Step.Time=2.5;
           SimulationSetup.Load.Step.Value=Rload;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=2.6;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.2;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=1.2;
           SimulationSetup.converter1.Pref.Values=[0.001;0.5];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,1.4];
           SimulationSetup.converter2.Pref.Values=[0,0.5];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[2.5-1e-3];

            
            SimulationSetup.timeVectorOperatingPointsEigenvalues=[2.499];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.05:0.3].';
            
  

    case 'SCR5_2MW_Rload_5percentStep'


        SCR=5;                   % short circuit ratio
        XR_ratio=10;             % X/R ratio  


        Scc=SCR*Sbase;           % short circuit capacity
        Xcct=Vbase^2/Scc;         % reactance 
        Rg_t=sqrt(Xcct^2/(XR_ratio^2+1));      
        Xg_t=Rg_t*XR_ratio;   

        % Thevenin equivalent grid 
        Zg_n2=complex(Rg_t,Xg_t);
        Rg=real(Zg_n2);         % Thevenin resistance
        Lg=imag(Zg_n2)/wg;
        Zgrid_pu=sqrt(Rg^2+Lg^2)/Zbase;


      SimulationSetup.Solver.timeSimulationEnd=3.5; 

           % load 
           SimulationSetup.Load.Base=1e2*Rload;
           SimulationSetup.Load.Step.Time=2.5;
           SimulationSetup.Load.Step.Value=1e2*Rload*1.05;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=2;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.2;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=1.2;
           SimulationSetup.converter1.Pref.Values=[0.00;0.2];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,1.4];
           SimulationSetup.converter2.Pref.Values=[0,0.0];

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[2.5-1e-3];

            
            SimulationSetup.timeVectorOperatingPointsEigenvalues=[2.499];
            SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.005:0.05].';


    case 'SCR5_200MW_Rload_5percentStep'
             SCR=5;                   % short circuit ratio
             XR_ratio=10;             % X/R ratio  
        
        
            Scc=SCR*Sbase;           % short circuit capacity
            Xcct=Vbase^2/Scc;         % reactance 
            Rg_t=sqrt(Xcct^2/(XR_ratio^2+1));      
            Xg_t=Rg_t*XR_ratio;   
    
            % Thevenin equivalent grid 
            Zg_n2=complex(Rg_t,Xg_t);
            Rg=real(Zg_n2);         % Thevenin resistance
            Lg=imag(Zg_n2)/wg;
            Zgrid_pu=sqrt(Rg^2+Lg^2)/Zbase;



           SimulationSetup.Solver.timeSimulationEnd=3; 

           % load 
           SimulationSetup.Load.Base=Rload*0.5;
           SimulationSetup.Load.Step.Time=2.5;
           SimulationSetup.Load.Step.Value=Rload*0.55;

           % grid disturbance 
           SimulationSetup.grid.voltageAmplitudeStep.Time=2;
           SimulationSetup.grid.voltageAmplitudeStep.Values=230e3*sqrt(2/3)*[1.01;0.92]; 

           SimulationSetup.grid.FrequencyStep.Time=10;
           SimulationSetup.grid.FrequencyStep.Values=[2*pi*50;2*pi*53];

           SimulationSetup.grid.PhaseStep.Time=4;
           SimulationSetup.grid.PhaseStep.Values=[0;(0.01*0.0174533)*1.06];
    
           % converter 1
           SimulationSetup.converter1.FreqCtrl.ActivationTime=0.2;
           SimulationSetup.converter1.VoltageCtrl.ActivationTime=0.8;

           SimulationSetup.converter1.Pref.Time=1.2;
           SimulationSetup.converter1.Pref.Values=[0.25;0.25];

           % converter 2
           SimulationSetup.converter2.Pref.Time=[0,0.2,0.21,0.22,0.3]; % [0,0.2]
           SimulationSetup.converter2.Pref.Values=[0,0.01,0.05,0.1,1]; % [0.5,1]

           % simulation and eigenvalue analysis points
           SimulationSetup.timeVectorOperatingPointsSimulation=[2.5-1e-3];

            
           SimulationSetup.timeVectorOperatingPointsEigenvalues=[2.499];
           SimulationSetup.timeVectorEigenvalueTrajectory=SimulationSetup.timeVectorOperatingPointsSimulation(1)+[0:0.005:0.05].';
            
   



    otherwise

        error('No case study was selected.')
end 



