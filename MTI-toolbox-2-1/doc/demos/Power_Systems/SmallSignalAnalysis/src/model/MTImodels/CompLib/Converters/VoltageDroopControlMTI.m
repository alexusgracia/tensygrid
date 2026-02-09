function [avrMTI] = VoltageDroopControlMTI(PrefixName,PrefixNumber,wf,Kdroop_V,Vpeak,Vqref)
%AVRMTI Summary of this function goes here
%   Detailed explanation goes here

% x1: xI_avr
% u1: Qref,  u2: Qmeas, 
% y1: Vdref, y2: Vqref

eqs={'0==-xp1-wf*x1+wf*u2';
     '0==-y1+Kdroop*(u1-x1)+Vpeak';
     '0==-y2+Vqref';    
     };

avrMTI=sym2dmss(eqs,0);

new=[wf,Kdroop_V,Vpeak,Vqref];
syms wf Kdroop Vpeak Vqref
old=[wf Kdroop Vpeak Vqref];

avrMTI=replaceSymbolicParameters(avrMTI,old,new);

avrMTI.stateName=[PrefixName+num2str(PrefixNumber)+"_xI_Vdroop"];
avrMTI.inputName=[PrefixName+num2str(PrefixNumber)+"_Qref",PrefixName+num2str(PrefixNumber)+"_Qmeas"];

avrMTI.algebraicName=[PrefixName+num2str(PrefixNumber)+"_Vdref",PrefixName+num2str(PrefixNumber)+"_Vqref"];
end

