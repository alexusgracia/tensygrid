function [ccLoop] = CurrentControlMTI(PrefixName,PrefixNumber,IdqRef,Vdq_pcc_local,Idq_local,omegaPLL,CC_Kp,CC_Ki,CC_Lf)
%CCLOOPMTI Creates a current control loop in dq-coordinates
%
% Input arguments: 
%       number- number of the converter 
%       Idq_ref -
%       Vdq_pcc_local - Vdq local of the PCC (string)
%       CC_Kp, CC_Ki - tuning parameter of the PI controller 
%       CC_Lf -     filter inductance
% Output arguments:
%       ccLoop - dmss-object of the current-controller 
%
% Future features: 
%

% u1-2: Idqref, u3-4: Vdq, u5-6: Idq, u7: omegaPLL

 % Current controller
  eqs={'xp1==Ki_CC*(u1-u5)';
       'xp2==Ki_CC*(u2-u6)';
       'y1==x1+Kp_CC*(u1-u5)-u6*Lf*u7+u3';
       'y2==x2+Kp_CC*(u2-u6)+u5*Lf*u7+u4';};

  ccLoop=sym2dmss(eqs,0);
  
      ccLoop.stateName=[PrefixName+num2str(PrefixNumber)+"_xId_CC",PrefixName+num2str(PrefixNumber)+"_xIq_CC"];
      ccLoop.inputName=[IdqRef,Vdq_pcc_local,Idq_local,omegaPLL];
      ccLoop.algebraicName=[PrefixName+num2str(PrefixNumber)+"_Ud",PrefixName+num2str(PrefixNumber)+"_Uq"];
       syms Ki_CC Kp_CC Lf w0
       
       
      ccLoop=replaceSymbolicParameters(ccLoop,[Ki_CC,Kp_CC,Lf,w0],[CC_Ki,CC_Kp,CC_Lf,2*pi*50]);

end

