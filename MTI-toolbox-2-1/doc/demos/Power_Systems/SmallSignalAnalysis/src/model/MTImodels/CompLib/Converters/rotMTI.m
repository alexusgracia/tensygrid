function [rot] = rotMTI(PrefixName,PrefixNumber,inRotName, trigConvAngle,Orientation)
%ROTMTI Rotation from DQ global into dq local 
% 

if size(inRotName)~=[1 2]
    error('inRotName must be of size 1 x 2.')
end
    % Rotation from global to local

if size(trigConvAngle)==[1 2] 
     % u1-2: dq input signal
     % u3: cos(delta_omega_local*t), u4: sin(delta_omega_local*t) 
     eqs={'y1==u3*u1+u4*u2';
              'y2==-u4*u1+u3*u2';};
elseif size(trigConvAngle)==[1 4]

    switch Orientation

        case 'dq'
            % u1-2: dq input signal
            % u3: cos(theta_local), u4: sin(theta_local) 
            % u5: cos(theta_global), u6: sin(theta_global)
        
            % Rotation from global to local
                   eqs={'y1==(u3*u5+u4*u6)*u1+(u4*u5-u3*u6)*u2';
                      'y2==-(u4*u5-u3*u6)*u1+(u3*u5+u4*u6)*u2';};

        case 'qd'

              % Rotation from global to local
                   eqs={'y1==(u3*u5+u4*u6)*u1-(u4*u5-u3*u6)*u2';
                      'y2==(u4*u5-u3*u6)*u1+(u3*u5+u4*u6)*u2';};

    end
else 
    error('trigConvAngle must be either of size 1 x 2 if the delta omega is used, or 1 x 4 if the angle difference between the local angle and global angle is taken.')
end

         rot=sym2dmss(eqs,0);
         rot.inputName=[inRotName,trigConvAngle];
         outName=[PrefixName+num2str(PrefixNumber)+"_"+inRotName];
         outName=strrep(outName,"D","d");
         outName=strrep(outName,"Q","q");
         rot.algebraicName=outName;
end

