
%% params
%% Base values
Sbase = 10e3;
Vbase = 380;
fbase = 50;
wbase = 2*pi*fbase;

nv = 9; % number of busess
num_inv = 3*3; 
Lines = 9;
w01 = 2*pi*50;
wn = ones(num_inv,1).*w01; % Nominal frequencies
Vn = 310.*ones(num_inv,1); % Nominal voltages

%% Single inverter parameter
% AC
Lf1 = 1.35e-3; %H 
nq1 = (3/100)*Vbase/Sbase;
%mp1 = (0.00002)*wbase/Sbase;
mp1 = (1/100)*wbase/Sbase;
Cf1 = 50e-6;   %F 
rf1 = 0.1; 
rc1 = 0.03;
%Lc1 = 0.35e-3; %H 
Lc1 = 0.0123;
% Voltage pi
Kpv1 = (1.25)/25;
Kiv1 = (10e3)/25;
% current pi
Kpc1 = (0.42)*(sqrt(3)*Vbase^2)/Sbase;
Kic1 = (640)*(2*wbase*Sbase/(sqrt(3)*Vbase^2));
F1 = 0.75;
fc1 = 5;
wc1 = 2*pi*fc1;
%DC
rdc1 = 28.86;
Cdc1 = 5;
Ldc1 = 1;
Vdc1 = 1000;


%% Matrix Values for all inverters

Lf = diag(Lf1.*ones(nv,1));
nq = diag(nq1.*ones(nv,1));
mp = diag(mp1.*ones(nv,1));
Cf = diag(Cf1.*ones(nv,1));
Kpv = diag(Kpv1.*ones(nv,1));
rf = diag(rf1.*ones(nv,1));
Kiv = diag(Kiv1.*ones(nv,1));
Lc = diag(Lc1.*ones(nv,1));
Kpc = diag(Kpc1.*ones(nv,1));
rc = diag(rc1.*ones(nv,1));
Kic = diag(Kic1.*ones(nv,1));
F = diag(F1.*ones(nv,1));
fc = diag(fc1.*ones(nv,1));
wc = diag(wc1.*ones(nv,1));
rdc = diag(rdc1.*ones(nv,1));
Cdc = diag(Cdc1.*ones(nv,1));
Ldc = diag(Ldc1.*ones(nv,1));
Vdc = Vdc1.*ones(nv,1);
Vtri = Vdc./2;
disp("Inverter parameters in workspace")














