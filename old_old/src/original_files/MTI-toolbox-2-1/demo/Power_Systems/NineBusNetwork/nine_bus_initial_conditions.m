

%% Network values in dq0
% Network data obtained from powerflow simulation of the nine-bus

% Nominal value for the sources
Vpeak = [310.1, 310.1, 310, 310.1, 310.1, 310, 310.1, 310.1 310];
Ipeak = [13.4 5.448 11.23 13.31 5.382 11.23 13.36 5.315 11.23];

% Phase angles
current_delta = (pi/180).*[-5.21 45.47 -13.82 -5.51 45.58 -13.83 -5.32 45.71 -13.82];
voltage_delta = (pi/180).*[0 0.1112 -0.2954 -0.0155 0.0964 -0.3037 -0.0155 0.0960 -0.2965]';

% Line Currents initial values in abc
fline0abc = [-2.262 1.621 -0.9812 -2.255 1.589 -1.013 -2.248 1.557 -1.046;
             1.932 -3.319 2.036 1.926 -3.258 2.097 1.920 -3.196 2.157;
             0.3305 1.698 -1.055 0.3289 1.669 -1.083 0.3279 1.64 -1.111];

% Initial Voltage at the nodes, apply Park Transfrom at t = 0.
V0bDQ = zeros(2,num_inv);
for i=1:num_inv
    aux = Vpeak(:,i)*(2/3).*[cos(0) cos(-2*pi/3) cos(2*pi/3); -sin(0) -sin(-2*pi/3) -sin(2*pi/3); 0.5 0.5 0.5]*[cos(voltage_delta(i));cos(voltage_delta(i)-2*pi/3);cos(voltage_delta(i)+2*pi/3)];
    V0bDQ(:,i) = aux(1:2);
end
V0bDQ = reshape(V0bDQ,[2*num_inv 1]);

% Initial Current injection at the nodes
I0DQ = zeros(2,num_inv);
for i=1:num_inv
    aux = Ipeak(:,i)*(2/3).*[cos(0) cos(-2*pi/3) cos(2*pi/3); -sin(0) -sin(-2*pi/3) -sin(2*pi/3); 0.5 0.5 0.5]*[cos(current_delta(i));cos(current_delta(i)-2*pi/3);cos(current_delta(i)+2*pi/3)];
    I0DQ(:,i) = aux(1:2);
end
I0DQ = reshape(I0DQ,[2*num_inv 1]);

% Line currents in DQ
fline0 = zeros(2,Lines);
for i=1:Lines
    fline0(:,i) = (2/3).*[cos(0) cos(-2*pi/3) cos(2*pi/3); -sin(0) -sin(-2*pi/3) -sin(2*pi/3)]*fline0abc(:,i);
end
fline0 = reshape(fline0,[2*Lines 1]);

%% Inverter values in dq0

% Set initial displacement as the voltage phase angle
delta0 = voltage_delta;

% Local io and vb initial values
V0bdq = zeros(size(V0bDQ));
I0dq = zeros(size(I0DQ));
for i=1:num_inv
    V0bdq(2*i-1:2*i,:) = [cos(voltage_delta(i)) -sin(voltage_delta(i));sin(voltage_delta(i)) cos(voltage_delta(i))]'*V0bDQ(2*i-1:2*i,:);
    I0dq(2*i-1:2*i,:) = [cos(voltage_delta(i)) -sin(voltage_delta(i));sin(voltage_delta(i)) cos(voltage_delta(i))]'*I0DQ(2*i-1:2*i,:);
end

% Clear voltage_delta, it is now in delta0
clear voltage_delta

% Take an initial guess on the inverter output voltage vo
V0dq = [310+round(rand(1,num_inv),1);zeros(1,num_inv)];
V0dq = reshape(V0dq,[2*num_inv 1]);

% Take an initial guess on the il current
diff_il0_io0 = zeros(2,num_inv);
diff_il0_io0(2,:) = -5.*ones(1,num_inv);
IL0dq = I0dq + reshape(diff_il0_io0,[2*num_inv 1]);

% Initial values for the instant power, filter initial conditions are set to zero
p0 = ([I0DQ(1:2:num_inv*2)].*V0bDQ(1:2:num_inv*2) + I0DQ(2:2:num_inv*2).*V0bDQ(2:2:num_inv*2));
q0 = (-[I0DQ(1:2:num_inv*2)].*V0bDQ(2:2:num_inv*2) + I0DQ(2:2:num_inv*2).*V0bDQ(1:2:num_inv*2));
wcom0 = -mp(1,1)*p0(1) + wn(1);

% DC circuit initialization idc0 and vdc0
Vdc = ones(num_inv,1).*Vdc(1);
idc0 = (3*310.*I0DQ(1:2:num_inv*2))./(Vdc); % 3phase-power over voltage
vdc0 = -rdc*idc0 + Vdc;
zdc0 = 1./vdc0;

% Inverter vl, v0ref and vlref
VL0dq = zeros(size(V0dq));
VL0dq(1:2:num_inv*2) = (rf*IL0dq(1:2:num_inv*2)) - Lf*wcom0*eye(num_inv)*IL0dq(2:2:num_inv*2) + V0dq(1:2:num_inv*2);
VL0dq(2:2:num_inv*2) = (rf*IL0dq(2:2:num_inv*2)) + Lf*wcom0*eye(num_inv)*IL0dq(1:2:num_inv*2) + V0dq(2:2:num_inv*2);

% Inital voref
voref0 = [kron(-nq*q0,[1 0]')] + [kron(Vn,[1 0]')];
dphi0 = voref0 - V0dq;

% Initial vlref0
vlref0 = VL0dq./kron(vdc0./(2*Vtri(1)),[1 1]');

% Initial ilref0
dgamma0 = zeros(2*num_inv,1);
ilref0 = dgamma0 - IL0dq;

% PI controller's inital phi and gamma
phi0 = kron(Kiv,eye(2))\(ilref0 - kron(Kpv,eye(2))*voref0 + kron(Kpv,eye(2))*V0dq - kron(F,eye(2))*I0dq - kron(diag(wn)*Cf,[0 -1;1 0])*V0dq);
gamma0 = kron(Kic,eye(2))\(vlref0 - kron(Kpc,eye(2))*ilref0 + kron(Kpc,eye(2))*IL0dq - kron(diag(wn)*Lf,[0 -1;1 0])*IL0dq);

%% In case multilinear with xtrigs and no constraints
for i=1:num_inv
     X0trig(2*i-1:2*i,:) = [cos(delta0(i)); sin(delta0(i))];
end

%% Set initial conditions estimates in the state vector

y0est = zeros(size(vars));
% Put initial states in their correct position in the state vector.
y0est(1:num_inv) = idc0;
y0est(num_inv+1:2*num_inv) = vdc0;
y0est(2*num_inv+1:3*num_inv) = delta0;
y0est(3*num_inv+1:4*num_inv) = zeros(size(num_inv,1)); % P0
y0est(4*num_inv+1:5*num_inv) = zeros(size(num_inv,1)); % Q0
y0est(5*num_inv+1:7*num_inv) = phi0;
y0est(7*num_inv+1:9*num_inv) = gamma0;
y0est(9*num_inv+1:11*num_inv) = IL0dq;
y0est(11*num_inv+1:13*num_inv) = V0dq;
y0est(13*num_inv+1:15*num_inv) = I0dq;
y0est(15*num_inv+1:15*num_inv+2*Lines) = fline0;
y0est(15*num_inv+2*Lines+1:15*num_inv+2*Lines+2*Nodes) = X0trig; 
y0est(15*num_inv+2*Lines+2*Nodes + 1 : 15*num_inv+2*Lines+2*Nodes+2*num_inv) = V0bDQ; 
y0est(15*num_inv+2*Lines+2*Nodes+2*num_inv+1:end) =zdc0;


disp("Initial conditions in workspace under y0est vector")