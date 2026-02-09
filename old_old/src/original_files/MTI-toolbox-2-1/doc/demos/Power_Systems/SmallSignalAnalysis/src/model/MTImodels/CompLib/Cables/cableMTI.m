function [cableSys] = cableMTI(PrefixName,PrefixNumber,cableModel,cableParam,coordSys,inputName)
%CABLEMTI - Creates an iMTI model of a desired cable section
%   
%    It is assumed that the current flows from the first node (first input voltage) to
%    the second node (second voltage).  
% 
%   Input arguments:
%       name - additional name if desired 
%       number 
%       coordSys - coordinate system, either 'abc' or 'dq'



% Features for the future:
%       ts       - discretization step size 
%       w0       - different angular frequency

ts=0; % currently only continuous-time systems

switch coordSys
    case 'abc'
        
         switch cableModel
            case 'RL'
                  eqs={'u1=R1*x1+L1*xp1+u4';
                       'u2=R1*x2+L1*xp2+u5';
                       'u3=R1*x3+L1*xp3+u6';
                      };
                
                cableSys=sym2dmss(eqs,ts);
                cableSys.stateName=[PrefixName+num2str(PrefixNumber)+"_Ia",PrefixName+num2str(PrefixNumber)+"_Ib",PrefixName+num2str(PrefixNumber)+"_Ic"];
                cableSys.inputName=inputName;

                syms R1 L1 
                old=[R1 L1];
             
                new=[cableParam];
                cableSys=replaceSymbolicParameters(cableSys,old,new);

             otherwise 
                 error('Only RL-models are currently implemented in abc-frame.')
         end

    case 'dq'
       
        switch cableModel
            case 'RL'

                % eqs={'u1=R*x1+L*xp1-w0*L*x2+u3';
                %        'u2=R*x2+L*xp2+w0*L*x1+u4'
                %       };

                eqs={'xp1=-R/L*x1+w0*x2+1/L*(u1-u3)';
                     'xp2=-w0*x1-R/L*x2+1/L*(u2-u4)';
                    };
                
                cableSys=sym2dmss(eqs,ts);
                cableSys.stateName=[PrefixName+num2str(PrefixNumber)+"_ID",PrefixName+num2str(PrefixNumber)+"_IQ"];
               
               % cableSys.inputName=["Vd"+num2str(number),"Vq"+num2str(number),"Vd_pcc"+num2str(number),"Vq_pcc"+num2str(number)];
               cableSys.inputName=inputName;
                %cableSys.algebraicName=["Vd_pcc"+num2str(number),"Vq_pcc"+num2str(number)]; % PCC voltage

                syms R L w0
                old=[R L w0];

                w0=2*pi*50;
                new=[cableParam, w0];
                cableSys=replaceSymbolicParameters(cableSys,old,new);
            case 'pi'
                error('pi model not implemented yet')
        end
end



end

