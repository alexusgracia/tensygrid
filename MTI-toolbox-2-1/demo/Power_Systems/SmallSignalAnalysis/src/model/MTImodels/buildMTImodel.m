%% Network model 
grid=gridMTI("grid","");

grid_RL=cableMTI("grid","",'RL',[Rg,Lg],'dq',["VD_g","VQ_g","load1_VD","load1_VQ"]);

switch SimulationSetup.model
    case 'GFM'
        load1=RloadMTI("load",1,'dq',[grid_RL.stateName,"converter1_ID","converter1_IQ"]);

    case 'GFM_GFL'
        load1=RloadMTI("load",1,'dq',[grid_RL.stateName,"converter1_ID","converter1_IQ","converter2_ID","converter2_IQ"]);
       
    case 'GFL'
        load1=RloadMTI("load",1,'dq',[grid_RL.stateName,"converter2_ID","converter2_IQ"]);    
end
%% MTI GFM model
converter1_Filter=cableMTI("converter",1,'RL',[Rconv1,Lconv1],"dq",["converter1_UD","converter1_UQ",load1.algebraicName]);
converter1_vsm=VirtualSynchronousMachineControlMTI("converter",1,H,K_damping,Sbase,wg);
converter1_Vdq=rotMTI("converter",1,load1.algebraicName,[converter1_vsm.stateName(2:3),grid.stateName(1:2)],'dq');
converter1_Idq=rotMTI("converter",1,converter1_Filter.stateName,[converter1_vsm.stateName(2:3),grid.stateName(1:2)],'dq');
converter1_PQ=PQcalculationMTI("converter",1,'dq',converter1_Vdq.algebraicName,converter1_Idq.algebraicName);
converter1_VoltageDroopControl=VoltageDroopControlMTI("converter",1,wf_PQ,Kdroop_V,Vpeak,0);
converter1_VA=virtualAdmittanceMTI("converter",1,Lv,Rv,wg,[converter1_VoltageDroopControl.algebraicName,converter1_Vdq.algebraicName]);
converter1_CC=CurrentControlMTI("converter",1,converter1_VA.algebraicName,converter1_Vdq.algebraicName,converter1_Idq.algebraicName,converter1_vsm.algebraicName(1),Kp_CC,Ki_CC,Lconv1);
converter1_Udq=INVrotMTI(converter1_CC.algebraicName,[converter1_vsm.stateName(2:3),grid.stateName(1:2)]);

% Connect 
 converter1=connect(converter1_Vdq,converter1_Filter,converter1_Idq,converter1_PQ,converter1_vsm,converter1_VA,converter1_CC,converter1_Udq,converter1_VoltageDroopControl);

%% MTI GFL model - current control


converter2_Filter=cableMTI("converter",2,'RL',[Rconv1,Lconv1],"dq",["converter2_UD","converter2_UQ",load1.algebraicName]);
converter2_PLL=PLLMTI("converter",2,'dq_theta',[load1.algebraicName,grid.stateName(1:2)],kp_PLL,ki_PLL,wg);
converter2_Vdq=rotMTI("converter",2,load1.algebraicName,[converter2_PLL.stateName(1:2),grid.stateName(1:2)],'dq');
converter2_Idq=rotMTI("converter",2,converter2_Filter.stateName,[converter2_PLL.stateName(1:2),grid.stateName(1:2)],'dq');
converter2_PQ=PQcalculationMTI("converter",2,'dq',load1.algebraicName,converter2_Filter.stateName);
converter2_PQctrl=PQcontrolMTI("converter",2,["converter2_Pref","converter2_Qref"],converter2_PQ.algebraicName,kp_P,ki_P,kp_Q,ki_Q);
converter2_CC=CurrentControlMTI("converter",2,converter2_PQctrl.algebraicName,converter2_Vdq.algebraicName,converter2_Idq.algebraicName,converter2_PLL.algebraicName(3),Kp_CC,Ki_CC,Lconv1);
converter2_Udq=INVrotMTI(converter2_CC.algebraicName,[converter2_PLL.stateName(1:2),grid.stateName(1:2)]);

converter2=connect(converter2_Filter,converter2_PLL,converter2_Vdq,converter2_Idq,converter2_PQ,converter2_PQctrl,converter2_CC,converter2_Udq);

%% Connect submodels 
switch SimulationSetup.model
    case 'GFL'
        MTImodel=connect(load1,grid_RL,grid,converter2);
    case 'GFM'
        MTImodel=connect(load1,grid_RL,grid,converter1);
    case 'GFM_GFL'
        MTImodel=connect(converter1,load1,grid_RL,grid);
        MTImodel=connect(MTImodel,converter2);  
end

pathMTImodels=which("buildMTImodel");