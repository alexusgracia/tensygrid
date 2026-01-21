function [convSys] = converterMTI(number,Vsync,PLL_Kp,PLL_Ki,inRotName)
%CONVERTERMTI Creates an MTI model of a converter
%
% Input arguments: 
%       number- number of the converter 
%       Vsync - voltage for synchronization (string)
%
% Output arguments:
%


if size(Vsync)~=[1 2]
    error('Vsync must be of size 1 x 2.')
end

ts=0;

% PLL
     eqs={'xp1==-(x3-Kp_PLL*(u1*x2+u2*x1))*y2'; %z1=cos(theta)
          'xp2==(x3-Kp_PLL*(u1*x2+u2*x1))*y1'; %z2=sin(theta)
          'xp3==-Ki_PLL*(u1*x2+u2*x1)';   %xI_PLL
          'x1==y1';
          'x2==y2';  %x4*x4+x5*x5-2*x5==0;
          'y3==(x3-Kp_PLL*(u1*x2+u2*x1)+2*pi*50)/(2*pi)'; % freq
          };
       PLL=sym2dmss(eqs,ts);
       PLL.stateName=["cos(theta_conv"+num2str(number)+")","sin(theta_conv"+num2str(number)+")","xI_PLL"+num2str(number)];
       PLL.inputName=Vsync;
       PLL.algebraicName=["cos(theta_conv"+num2str(number)+")_algb","sin(theta_conv"+num2str(number)+")_algb","freq_conv"+num2str(number)];
       syms Kp_PLL Ki_PLL
       old=[Kp_PLL,Ki_PLL];
       new=[PLL_Kp,PLL_Ki];
       PLL=replaceSymbolicParameters(PLL,old,new);


         % Rotation from global to local
         eqs={'y1==u3*u1-u4*u2';
              'y1==u4*u1+u3*u2';}
         rot=sym2dmss(eqs,ts);
         rot.inputName=[inRotName,trigConvAngle];
         rot.algebraicName=inRotName'+["_local";"_local"];

    
      'y3==x3*y1-x4*y2'; % Vd local
      'y4==x4*y1+x3*y2'; % Vq local 
      %'y6==Kp_PLL*1/3*((-2*y1+y2+y3)*x5+sqrt(3)*(y2-y3)*x4)+w0+x6'; % w
      'y5==x3*x1-x4*x2'; % Id local
      'y6==x4*x1+x3*x2';  % Iq local 
       % Current controller
      'xp6==Ki_CC*(u3-y5)';
      'xp7==Ki_CC*(u4-y6)';
      'y7==x6+Kp_CC*(u3-y5)-y6*Lf*w0+y3';
      'y8==x7+Kp_CC*(u4-y6)+y5*Lf*w0+y4';
      %Inverse Park transformation   
       'y9==x3*y7+x4*y8';
       'y10==-x4*y7+x3*y8';
       'y13==(x5-Kp_PLL*(y1*x4+y2*x3)+2*pi*50)/(2*pi)';
        





outputArg1 = inputArg1;
conSys = inputArg2;
end

