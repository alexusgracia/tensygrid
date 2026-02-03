function [Rload] = RloadMTI(PrefixName,PrefixNumber,coorSys,inputName)
%RLOADMTI A resistive load connected to ground. 
%
% Input arguments:
%       inputName: string of the inputName
%       outputName: string of the output signal name, e.g. "V" or "I"
%       outputNameSubscript: string of the subscript of the output signal name, e.g. "conv"
%       outputNameNumber: number of the unit
%
% Output arguments:
%       Rload: dmss object 
%
% Date: 27.01.2025
% Author: Christoph Kaufmann

% if length(inputName) ==3
%     coorSys='abc';
% elseif length(inputName)==2 
%     coorSys='dq';
% else 
%     error('inputName must be a string of length of either 3 for abc, or 2 for dq-frame.')
% end 

switch coorSys
    % case 'abc'
        % Gs=eye(3); % Voltage is the input, current the output
        % Gphi=1/Rparam*eye(3); 
        % 
        % Rload=dmss(mss(CPN1(Gs,[]),CPN1(Gs,Gphi)));
        % Rload.algebraicName=['Ia_load_'+outputNameSubscript+num2str(outputNameNumber),'Ib_load_'+outputNameSubscript+num2str(outputNameNumber),'Ic_load_'+outputNameSubscript+num2str(outputNameNumber)]

    case 'dq'
        
        numberInputs=length(inputName)/2;
        Gs=[eye(length(inputName)),ones(length(inputName),1)].';
        D=[];

        for k=1:numberInputs
            D=[D,eye(2)];            
        end

  
        Gphi=D; 
        Rload=dmss(mss(CPN1(Gs,[]),CPN1(Gs,Gphi)));
        Rload.algebraicName=[PrefixName+num2str(PrefixNumber)+'_VD',PrefixName+num2str(PrefixNumber)+'_VQ'];

end
    
Rload.inputName=[inputName,PrefixName+num2str(PrefixNumber)+'_Rvariable'];
end

