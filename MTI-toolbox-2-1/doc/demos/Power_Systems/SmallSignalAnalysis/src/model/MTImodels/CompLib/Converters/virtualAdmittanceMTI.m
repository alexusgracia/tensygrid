function [VAmti] = virtualAdmittanceMTIv2(PrefixName,PrefixNumber,Lv,Rv,w0,inputName)
%VIRTUALADMITTANCEMTI Summary of this function goes here
%   Detailed explanation goes here

%u1-2: Vdqref, u3-4: Vdq

if size(inputName)~=[1 4]
    error('inputName must be of size 1 x 4, containing Vdqref and Vdq.')
end

eqs={'xp1==-Rv/Lv*x1+w0*x2+(u1-u3)';
     'xp2==-Rv/Lv*x2-w0*x1+(u2-u4)';
     'y1==1/Lv*x1'; 
     'y2==1/Lv*x2';
     };

VAmti=sym2dmss(eqs,0);
new=[Lv,Rv,w0];
syms Lv Rv w0
old=[Lv,Rv,w0];

VAmti=replaceSymbolicParameters(VAmti,old,new);

VAmti.inputName=inputName;
VAmti.stateName=[PrefixName+num2str(PrefixNumber)+"_xId_VA",PrefixName+num2str(PrefixNumber)+"_xIq_VA"];
VAmti.algebraicName=[PrefixName+num2str(PrefixNumber)+"_Idref",PrefixName+num2str(PrefixNumber)+"_Iqref"];
end

