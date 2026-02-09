function [vsmMTI] = VirtualSynchronousMachineControlMTI(PrefixName,PrefixNumber,H,K_damping,Sbase,w0)
%VSMMTI Summary of this function goes here
%   Detailed explanation goes here

% u1: Pref, u2: Pmeas
% u3: deblockingSignal 

% y: omegaVSM

eqs={  '0==-xp1-K_damping/(2*H)*x1+1/(2*H)*u1-1/(2*H*Sbase)*u2+K_damping/(2*H)';
       '0==-xp2-w0*(x1)*x3'; % cos(theta)
      '0==-xp3+w0*(x1)*x2'; % sin(theta)
      '0==-xp4+w0*(x1)'; % theta
      '0==-y1+w0*(x1)'; % omega
    };

   % 'xp1==-K_damping/(2*H)*x1+1/(2*H)*u1-1/(2*H*Sbase)*u2+Kdamping/(2*H)';
   % 'xp2==-w0*x1*x3*u3+w0*x3*u3-w0*x3'; % cos(theta)
   % 'xp3==w0*x1*x2*u3-w0*x2*u3+w0*x2'; % sin(theta)
   % 'xp4==+w0*x1*u3-w0*u3+w0'; % theta

Fs=[1 0 0 0 1 0 0 1 0 0 1 0; %x1
    0 0 0 0 0 0 0 1 1 1 0 0; %x2   
    0 0 0 0 1 1 1 0 0 0 0 0; %x3
    0 0 0 0 0 0 0 0 0 0 0 0; %x4
    0 1 0 0 0 0 0 0 0 0 0 0; %u1
    0 0 1 0 0 0 0 0 0 0 0 0; %u2
   % 0 0 0 0 1 1 0 1 1 0 1 1; %u3
];

K=K_damping; Sb=Sbase; w=w0;


Fphi=[-K/(2*H), 1/(2*H), -1/(2*H*Sb),K/(2*H), 0,0,0,0,0,0,0,0;
        0     , 0     ,   0,    0    ,-w,w,-w,0,0,0,0,0; 
       0     , 0     ,   0,    0    ,0,0,0,w,-w,w,0,0 ;
        0     , 0     ,   0,    w    ,0,0,0,0,0,0,w,-w  ;];

Gphi=[   0     , 0     ,   0,    w    ,0,0,0,0,0,0,w,-w  ;];

vsmMTI=sym2dmss(eqs,0);%dmss(mss(CPN1(Fs,Fphi),CPN1(Fs,Gphi)));


% vsmMTI=sym2dmss(eqs,0);
% 
new=[K_damping,H,w0,Sbase];
syms K_damping H w0 Sbase
old=[K_damping,H,w0,Sbase];
% 
vsmMTI=replaceSymbolicParameters(vsmMTI,old,new);

vsmMTI.inputName=[PrefixName+num2str(PrefixNumber)+"_Pref",PrefixName+num2str(PrefixNumber)+"_Pmeas"];
vsmMTI.stateName=[PrefixName+num2str(PrefixNumber)+"_xI_vsm",PrefixName+num2str(PrefixNumber)+"_cos(thetaVSM)",PrefixName+num2str(PrefixNumber)+"_sin(thetaVSM)",PrefixName+num2str(PrefixNumber)+"_thetaVSM"];
vsmMTI.algebraicName=[PrefixName+num2str(PrefixNumber)+"_omegaVSM"];
end

