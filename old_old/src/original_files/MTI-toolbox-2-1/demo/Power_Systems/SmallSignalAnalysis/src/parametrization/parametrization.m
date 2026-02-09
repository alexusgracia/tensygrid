%% Network parameters

% Base values 
Sbase=100e6;              % base apparent power (VA)
Vbase=230e3;              % nominal grid voltage (V)
Vpeak=Vbase*sqrt(2/3);
Ibase=Sbase/Vbase;
fg=50;                    % nominal grid frequency (Hz)
freq=fg;
wg=2*pi*fg;               % nom. angular grid freq. (rad/s)


% Grid parameters
Zbase=Vbase^2/Sbase;      % Base impedance  
XR_ratio=10;              % X/R ratio  

SCR2=5;                   % short circuit ratio
Scc2=SCR2*Sbase;           % short circuit capacity

Xcct=Vbase^2/Scc2;         % reactance 
Rg_t=sqrt(Xcct^2/(XR_ratio^2+1));      
Xg_t=Rg_t*XR_ratio;             

% Thevenin equivalent grid 
Zg_n2=complex(Rg_t,Xg_t);
Rg=real(Zg_n2);         % Thevenin resistance
Lg=imag(Zg_n2)/wg;
Zgrid_pu=sqrt(Rg^2+Lg^2)/Zbase;

Vdc_ref=400e3; 

Vpeak=Vbase/sqrt(3)*sqrt(2);     % Peak voltage
Xn=Vbase^2/Sbase;                % Base impedance
Ln=Xn/wg;                        % Base inductance
Ipeak = Sbase/Vbase*sqrt(2/3);

R_f_c=0.002*Xn;                   % Converter grid coupling filter resistance
X_f_c=0.07; 

L_f_c=X_f_c*Zbase/wg;
Z_filt=sqrt(R_f_c^2+X_f_c^2);
Zfilt_pu=Z_filt/Xn;
XR_filt=X_f_c/R_f_c;

%% Grid-forming parameters
% VSM 
H=1;
K_damping=50;
% VSM 2

%  % % VA 1
% Lv=L_f_c
% Rv=R_f_c

% VA 
Kva=0.02;

Lv=Kva*Zbase/wg;
tau_va=2*1e-3*56;
Rv=Lv/tau_va;
z_va=sqrt(Rv^2+(wg*Lv)^2)/Zbase;
XR_va=Lv*wg/Rv;

% Q droop gfor
wf_PQ=20;
tau_droop_u=1/wf_PQ;
mq_droop=0.2;
qddd=1/mq_droop;
Kdroop_V=mq_droop*Vpeak/Sbase;

%
% PLL
    Ki_PLL=10;
    Kp_PLL=1;

% Filter 
Rconv1=R_f_c;
Lconv1=L_f_c;


% Current controller 
  fsw=1e3;        
  tau_c=1/fsw; % time constant of the current control loop 10 times faster
  Kp_CC=Lconv1/tau_c;
  Ki_CC=Rconv1/tau_c;

% Load resistance
Pload=2*Sbase/3;
Rload=Vbase^2/Pload;


% Reference values GFM
Pref=1.375e6/Sbase;
Qref=Sbase*0.001;

%% Grid-following tuning

% pq control
tau_p=0.1;
kp_P=1e-3*tau_c/tau_p;
ki_P=1e-3*1/tau_p;
tau_q=0.1;
kp_Q=1e-3*tau_c/tau_q;
ki_Q=1e-3*1/tau_q;

% pll

ts_pll_n1=0.4;
xi_pll_n1=0.707;
omega_pll_n1=4/(ts_pll_n1*xi_pll_n1);
kp_PLL=xi_pll_n1*2*omega_pll_n1/Vpeak;
tau_pll_n1=2*xi_pll_n1/omega_pll_n1;
ki_PLL=kp_PLL/tau_pll_n1;

% w-P droop
tau_droop_f=1/20;
m_p_droop=20;
k_droop_f=m_p_droop*Sbase/wg;

% Q-V droop
tau_droop_u=1/20;
mq_droop=20;
k_droop_u_gfol=mq_droop*Sbase/Vpeak;
