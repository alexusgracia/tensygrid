function [PQcalc] = PQcalculationMTI(PrefixName,PrefixNumber,coorSys,VName,IName)
%PQCALCMTI Creates an MTI object of the PQ calculation  
%   
% Input arguments:
%       coorSys: choose between "abc" or "dq" 
%       VName: string of the inputName related to V
%       IName: string of the inputName related to I
%       OutputNamePrefix: string of the prefix for the OutputName
%       OutputNameNumber: number of the unit
%
% Output arguments:
%       PQcalc: dmss object 
%
% Future features 
%       - implementation for dq0
%
% Date: 27.01.2025
% Author: Christoph Kaufmann

if not(isstring(VName)) || not(isstring(IName))
    error ('VName and IName must be strings.')
end

switch coorSys 

    case 'abc'
          if any(size(VName)~=[1 3]) || any(size(IName)~=[1 3]) 
            error('VName and IName must be of size 1 x 3.')
           end

        Gs=[1 0 0 0 0 0 1 1 0; % u1=Va
            0 1 0 1 0 0 0 0 1; % 
            0 0 1 0 1 1 0 0 0; % u3=Vc
            1 0 0 1 1 0 0 0 0; % u4=Ia
            0 1 0 0 0 1 1 0 0; % 
            0 0 1 0 0 0 0 1 1; % u5=Ic
            ] ;

        Gphi=[1;1/sqrt(3)].*[1 1 1 0 0 0 0 0 0;
                             0 0 0 1 -1 1 -1 1 -1];

        PQcalc=dmss(mss(CPN1(Gs,[]),CPN1(Gs,Gphi)));

    case 'dq'   
        if any(size(VName)~=[1 2]) || any(size(IName)~=[1 2]) 
            error('VName and IName must be of size 1 x 2.')
        end
        
         % Fs=[0,0,0,0; %x1 
         %     0,0,0,0; %x2 
         %     1, 0,0,1; % Vd
         %     0, 1,1,0;  % Vq
         %     1, 0,1,0;  % Id
         %     0  1,0,1] % Iq 
         % Fphi=3/2*[1, 1, 0, 0;
         %           0, 0, 1, -1];

        ts=0;
        
        % assumes that there is no 0-component!!
        eqs={'y1=3/2*(u1*u3+u2*u4)';
             'y2=3/2*(u2*u3-u1*u4)'};
        PQcalc=sym2dmss(eqs,ts);
            

    case 'dq0'
         % if any(size(VName)~=[1 2]) || any(size(IName)~=[1 2]) 
         %  error('VName and IName must be of size 1 x 3.')
         % end

         error('dq0 is not implemented yet ')

end 

PQcalc.inputName=[VName,IName];
PQcalc.algebraicName=[PrefixName+num2str(PrefixNumber)+"_Pmeas",PrefixName+num2str(PrefixNumber)+"_Qmeas"];
end


