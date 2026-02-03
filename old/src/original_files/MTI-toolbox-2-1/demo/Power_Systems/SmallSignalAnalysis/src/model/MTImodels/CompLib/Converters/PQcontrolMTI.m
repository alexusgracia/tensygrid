function [PQctrl] = PQcontrolMTI(PrefixName,PrefixNumber,PQreferenceName,PQcalculationName,Kp_P,Ki_P,Kp_Q,Ki_Q)
%PQcontMTI Creates a MTI object of a PQ-controller with Idqref as outputs
%   
% Input parameters:
%   PrefixName: prefix of the variable name
%   PrefixName: number of the device unit associated with the prefix of the variable name
%   Ki_P,Kp_P,Ki_Q,Ki_Q: controller gains
% Output parameters:
%
% Future features: 
%    - choose between mss and dmss object 
% KauChr, 10.06.2025


if size(PQreferenceName)~=[1 2]
    error('PQrefName must be of size 1 x 2.')
end

if size(PQcalculationName)~=[1 2]
    error('PQcalcName must be of size 1 x 2.')
end

% state equations
Fs=[0 0 0 0;% xId
    0 0 0 0;% xIq
    1 0 0 0;% Pref
    0 1 0 0;% Pmeas
    0 0 1 0;% Qref
    0 0 0 1;];% Qmeas
Fphi=[Ki_P, -Ki_P, 0, 0;       %xpId
       0,     0, Ki_Q,-Ki_Q];  %xpIq
% output equations
Gs=[0 0 1 0 0 0;% xId
    0 0 0 0 0 1;% xIq
    1 0 0 0 0 0;% Pref
    0 1 0 0 0 0;% Pmeas
    0 0 0 1 0 0;% Qref
    0 0 0 0 1 0;];% Qmeas

Gphi=[Kp_P, -Kp_P, 1, 0, 0, 0;       %Idref
       0,     0,  0, -Kp_Q,Kp_Q, -1]; %Iqref

PQctrl=mss(CPN1(Fs,Fphi),CPN1(Gs,Gphi));

PQctrl=dmss(PQctrl); 

PQctrl.stateName=[PrefixName+num2str(PrefixNumber)+"_xI_Pctrl",PrefixName+num2str(PrefixNumber)+"_xI_Qctrl"];
PQctrl.inputName=[PQreferenceName(1),PQcalculationName(1),PQreferenceName(2),PQcalculationName(2)];
PQctrl.algebraicName=[PrefixName+num2str(PrefixNumber)+"_Idref",PrefixName+num2str(PrefixNumber)+"_Iqref"];

end

