function [PLL] = PLLMTI(PrefixName,PrefixNumber,coorSys,VsyncName,PLL_Kp,PLL_Ki,wnom)
%PLLMTI Creates an iMTI model of a PLL
%
% Input arguments: 
%       number- number of the converter 
%       Vsync - voltage for synchronization (string)
%       wnom - nominal angular frequency, e.g. 2*pi*50
% Output arguments:
%
%
% Version log:
% v3 - removed negative sign in vq input signal

if nargin<3 
    fprintf('No gain values and nominal grid frequency provided for the PLL, assuming Kp=1, Ki=10, and 50 Hz.')
    PLL_Kp=1;
    PLL_Ki=10;
    wnom=2*pi*50;
elseif nargin<5
    fprintf('No nominal grid frequency provided for the PLL, assuming 50 Hz.')
    wnom=2*pi*50;
else
    
end 

ts=0;


% PLL
switch coorSys
    case 'abc'
        if size(VsyncName)~=[1 3] 
            error('Vsync must be of size 1 x 3.')
         end
    % q-axis is aligned with the a-phase, meaning Vq is driven to zero and Vd is the amplitude of the voltage    
    eqs={'0==-xp1+(x3+Kp_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*(x2))+w0)*(-y2)'; % x1=z1=cos(theta)
     '0==-xp2+(x3+Kp_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*(x2))+w0)*y1'; % x2=z2=sin(theta) 
     '0==-xp3+(Ki_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*(x2)))';   % x3=z3=xI_PLL
     '0==-xp4+(Kp_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*x2)+w0+x3)';
     '0==y1-x1';
     '0==y2-x2'; %x1*y1+x2*y2-1==0;
     'y3==(Kp_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*x2)+w0+x3)'; % % anngular freq    
     }; 

    case 'dq_delta' % uses the angle difference \delta=\tehta_i-\theta_{grid}

        if size(VsyncName)~=[1 2] 
            error('Vsync must be of size 1 x 2.')
         end



    eqs={'xp1==-(x3+Kp_PLL*(-u1*x2+u2*x1))*y2'; %z1=cos(theta)
          'xp2==(x3+Kp_PLL*(-u1*x2+u2*x1))*y1'; %z2=sin(theta)
          'xp3==+Ki_PLL*(-u1*x2+u2*x1)';   %xI_PLL
          'x1==y1';
          'x2==y2';  %x4*x4+x5*x5-2*x5==0;
          'y3==(x3+Kp_PLL*(-u1*x2+u2*x1)+w0)'; % anngular freq
          };

    case 'dq_theta' % uses \theta_i and theta_grid

        if size(VsyncName)~=[1 4] 
            error('Vsync must be of size 1 x 4, i.e., [VD, VQ, cos(theta_grid), sin(theta_grid)].')
         end

        eqs={'xp1==-(x3+Kp_PLL*((u4*u1+u3*u2)*x1+(-u3*u1+u4*u2)*x2)+w0)*y2'; %z1=cos(theta)
          'xp2==(x3+Kp_PLL*((u4*u1+u3*u2)*x1+(-u3*u1+u4*u2)*x2)+w0)*y1'; %z2=sin(theta)
          'xp3==Ki_PLL*((u4*u1+u3*u2)*x1+(-u3*u1+u4*u2)*x2)';   %xI_PLL
          '0==-xp4+(x3+Kp_PLL*((u4*u1+u3*u2)*x1+(-u3*u1+u4*u2)*x2)+w0)'; %thetaPLL
          'x1==y1';
          'x2==y2';  %x4*x4+x5*x5-2*x5==0;
          'y3==(x3+Kp_PLL*((u4*u1+u3*u2)*x1+(-u3*u1+u4*u2)*x2)+w0)'; % anngular freq
          };
%'0==-xp4+(Kp_PLL*1/3*((-2*u1+u2+u3)*x1+sqrt(3)*(u2-u3)*x2)+w0+x3)'; %thetaPLL
        



    case 'qd'
        
            if size(VsyncName)~=[1 2] 
            error('Vsync must be of size 1 x 2.')
         end



         eqs={'xp1==-(x3-Kp_PLL*((u3*u1-u4*u2)*x1+(u3*u2+u4*u1)*x2)+w0)*y2'; %z1=cos(theta)
          'xp2==(x3-Kp_PLL*((u3*u1-u4*u2)*x1+(u3*u2+u4*u1)*x2)+w0)*y1'; %z2=sin(theta)
          'xp3==-Ki_PLL*((u3*u2+u4*u1)*x2+(u3*u1-u4*u2)*x1)';   %xI_PLL
          'x1==y1';
          'x2==y2';  %x4*x4+x5*x5-2*x5==0;
          'y3==(x3-Kp_PLL*((u3*u1-u4*u2)*x1+(u3*u2+u4*u1)*x2)+w0)'; % anngular freq
          };



end 
       PLL=sym2dmss(eqs,ts);
       PLL.stateName=[PrefixName+num2str(PrefixNumber)+"_cos(thetaPLL)",PrefixName+num2str(PrefixNumber)+"_sin(thetaPLL)",PrefixName+num2str(PrefixNumber)+"_xI_PLL",PrefixName+num2str(PrefixNumber)+"_thetaPLL"];
       PLL.inputName=VsyncName;
       PLL.algebraicName=[PrefixName+num2str(PrefixNumber)+"_cos(thetaPLL)_algb",PrefixName+num2str(PrefixNumber)+"_sin(thetaPLL)_algb",PrefixName+num2str(PrefixNumber)+"_omegaPLL"];
       syms Kp_PLL Ki_PLL w0
       old=[Kp_PLL,Ki_PLL,w0];
       new=[PLL_Kp,PLL_Ki,wnom];
       PLL=replaceSymbolicParameters(PLL,old,new);
end

