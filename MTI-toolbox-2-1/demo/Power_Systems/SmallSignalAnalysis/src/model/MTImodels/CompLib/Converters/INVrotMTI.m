function [rot] = INVrotMTI(inRotName, trigConvAngle)
%INVROTMTI Inverse rotation from dq local to DQ global 
% 
% Improvements: better hints for the sequence of inputNames?. 
% 


if size(inRotName)~=[1 2]
    error('inRotName must be of size 1 x 2.')
end

if size(trigConvAngle)~=[1 2]
    error('trigConvAngle must be of size 1 x 2.')
end



 % inverse rotation from global to local
if size(trigConvAngle)==[1 2] 
     % u1-2: dq input signal
     % u3: cos(delta_omega_local*t), u4: sin(delta_omega_local*t) 
     eqs={'y1==u3*u1+u4*u2';
          'y2==-u4*u1+u3*u2';};
elseif size(trigConvAngle)==[1 4]
    % u1-2: dq input signal
    % u3: cos(theta_local), u4: sin(theta_local) 
    % u5: cos(theta_global), u6: sin(theta_global)

          eqs={'y1==(u3*u5+u4*u6)*u1-(u4*u5-u3*u6)*u2';
              'y2==(u4*u5-u3*u6)*u1+(u3*u5+u4*u6)*u2';};
else 
    error('trigConvAngle must be either of size 1 x 2 if the delta omega is used, or 1 x 4 if the angle difference between the local angle and global angle is taken.')
end


         rot=sym2dmss(eqs,0);
         rot.inputName=[inRotName,trigConvAngle];
         outName=inRotName;
         outName=strrep(outName,"d","D");
         outName=strrep(outName,"q","Q");
         rot.algebraicName=outName;
         %varName=extract(rot.algebraicName(1),1);
         %rot.algebraicName=rot.algebraicName.replace("_local","");
         %rot.algebraicName=rot.algebraicName.replace(varName,"D");
end

