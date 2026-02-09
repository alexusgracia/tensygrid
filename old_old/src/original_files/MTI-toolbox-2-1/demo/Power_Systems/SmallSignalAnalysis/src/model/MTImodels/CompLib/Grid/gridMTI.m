function [gridMTI] = gridMTIv2(PrefixName,PrefixNumber)
%GRIDMTI Summary of this function goes here
%   


eqs={'0==-xp1-u1*x2';
     '0==-xp2+u1*x1';
     '0==-xp3+u1'};

gridMTI=sym2dmss(eqs,0);

gridMTI.stateName=[PrefixName+num2str(PrefixNumber)+"_cos(thetaGRID)",PrefixName+num2str(PrefixNumber)+"_sin(thetaGRID)",PrefixName+num2str(PrefixNumber)+"_thetaGRID"];
gridMTI.inputName=[PrefixName+num2str(PrefixNumber)+"_omega"];
end

