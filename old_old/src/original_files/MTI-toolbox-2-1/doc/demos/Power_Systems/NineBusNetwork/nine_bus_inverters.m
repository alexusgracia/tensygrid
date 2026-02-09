%% Power Controller variables
vo = symDQvars('vo',num_inv);
io = symDQvars('io',num_inv);
il = symDQvars('il',num_inv);
omega = symSingleVars('omega',num_inv);
Q = symSingleVars('Q',num_inv);
delta = symSingleVars('delta',num_inv);
voref = symDQvars('voref',num_inv);     % vo voltage reference

%% Voltage Controller variables
phi = symDQvars('phi',num_inv);
ilref = symDQvars('ilref',num_inv);

%% Current Controller variables
gamma = symDQvars('gamma',num_inv);
vlref = symDQvars('vlref',num_inv);

%% Swtiching Mode variables
vl = symDQvars('vl',num_inv);

%% Bus voltages in dq
vbdq = symDQvars('vb',num_inv);

%% State Equations
syms wcom(t)
% DC - multilinear already in the demo
eq2_nl = diff(vdc,t) == Cdc\(idc-P./vdc); % nonlinear version
% Power controller 
eq4 = diff(delta,t)== -mp*P + wn-wcom;
   p = [io(1:2:num_inv*2)].*vo(1:2:num_inv*2) + io(2:2:num_inv*2).*vo(2:2:num_inv*2); 
   q = -[io(1:2:num_inv*2)].*vo(2:2:num_inv*2) + io(2:2:num_inv*2).*vo(1:2:num_inv*2);
eq5 = diff(P,t)== -wc*P + wc*p;  
eq6 = diff(Q,t)== -wc*Q + wc*q; 
% Voltage Controller
eq7 = diff(phi,t)== voref - vo; 
% Current Controller 
eq8 = diff(gamma,t) == ilref - il; 
% Switching mode
eq9 = vl == (vlref.*kron(vdc,ones(2,1)))./kron(2*Vtri',[1 1])'; 

% LCL filter, differential equations for (il,vo,io) in dq.
eq10 = diff(il(1:2:num_inv*2),t) == -Lf\(rf*il(1:2:num_inv*2)) + wcom*eye(num_inv)*il(2:2:num_inv*2) + Lf\(vl(1:2:num_inv*2)-vo(1:2:num_inv*2));
eq11 = diff(il(2:2:num_inv*2),t) == -Lf\(rf*il(2:2:num_inv*2)) - wcom*eye(num_inv)*il(1:2:num_inv*2) + Lf\(vl(2:2:num_inv*2)-vo(2:2:num_inv*2));
eq12 = diff(vo(1:2:num_inv*2),t) == wcom*eye(num_inv)*vo(2:2:num_inv*2) + Cf\(il(1:2:num_inv*2)-io(1:2:num_inv*2));
eq13 = diff(vo(2:2:num_inv*2),t) == -wcom*eye(num_inv)*vo(1:2:num_inv*2) + Cf\(il(2:2:num_inv*2)-io(2:2:num_inv*2));
eq14 = diff(io(1:2:num_inv*2),t) == -Lc\(rc*io(1:2:num_inv*2)) + wcom*eye(num_inv)*io(2:2:num_inv*2) + Lc\(vo(1:2:num_inv*2)-vbdq(1:2:num_inv*2));
eq15 = diff(io(2:2:num_inv*2),t) == -Lc\(rc*io(2:2:num_inv*2)) - wcom*eye(num_inv)*io(1:2:num_inv*2) + Lc\(vo(2:2:num_inv*2)-vbdq(2:2:num_inv*2));


%% Output equations for the controller blocks
% Power Controller
oeq1 = voref == [kron(-nq*Q,[1 0]')] + [kron(Vn,[1 0]')];
% Voltage Controller 
oeq2 = ilref == kron(Kiv,eye(2))*phi + kron(Kpv,eye(2))*voref - kron(Kpv,eye(2))*vo + kron(F,eye(2))*io + kron(diag(wn)*Cf,[0 -1;1 0])*vo;
% Current Controller 
oeq3 = vlref == kron(Kic,eye(2))*gamma + kron(Kpc,eye(2))*ilref - kron(Kpc,eye(2))*il + kron(diag(wn)*Lf,[0 -1;1 0])*il;

%% Substitutions to eleminate all reference values
% An aditional s is a added to every equation name where a substitution was made.

% Reference voltage and current
eq7s = subs(eq7,voref,rhs(oeq1));    % phi -> voref
oeq2s = subs(oeq2,voref,rhs(oeq1));  % ilref -> voref
eq8s = subs(eq8,ilref,rhs(oeq2s));   % gamma -> ilref
oeq3s = subs(oeq3,ilref,rhs(oeq2s)); % vlref -> ilref
eq9s = subs(eq9,vlref,rhs(oeq3s));   % vl -> vlref


%% Local dqi to global DQi transformation
% There are nonlinear, multilinear and linearized transformations.

% Bus voltages in DQ variables, same for all cases.
vbDQ = symDQvars('vbDQ',num_inv);

% Nonlinear rotation matrix
TDQ_nl=[];
for i=1:num_inv
    T = [cos(delta(i)),-sin(delta(i));sin(delta(i)) cos(delta(i))]; % fDQ = T*fdq
    TDQ_nl = blkdiag(TDQ_nl,transpose(T));
end

% These equations will have the trigonometric functions.
eq14s_nl = subs(eq14,vbdq,TDQ_nl*vbDQ);
eq15s_nl = subs(eq15,vbdq,TDQ_nl*vbDQ);

% Linearized, look into equations (37) to (41) in the Nagaraju Pogaku - Modeling, Analysis and Testing of Autonomous
% Operation of an Inverter-Based Microgrid.
% 
% Linearized rotation matrix
% Ts=[];
% for i=1:num_inv
%     T = [cos(delta0(i)),-sin(delta0(i));sin(delta0(i)), cos(delta0(i))]; 
%     Ts = blkdiag(Ts,T);
% end
% Tc=[];
% for i=1:num_inv
%     T = [-sin(delta0(i)),-cos(delta0(i));cos(delta0(i)) -sin(delta0(i))]*I0dq(2*i-1:2*i,:); 
%     Tc = blkdiag(Tc,T);
% end
% Tv=[];
% for i=1:num_inv
%     T = [-sin(delta0(i)),cos(delta0(i));-cos(delta0(i)) -sin(delta0(i))]*V0bDQ(2*i-1:2*i,:); 
%     Tv = blkdiag(Tv,T);
% end
% 
% % Linearization of vDQ = T*vdq.
% aux_eq = vbdq == V0bdq + Tv*(delta-delta0) + Ts'*(vbDQ - V0bDQ);
% eq14s = subs(eq14,vbdq,rhs(aux_eq));
% eq15s = subs(eq15,vbdq,rhs(aux_eq));

% Multilinear with no constraints.
% Introducing dummy varaibles as,
%       z1 = cos(delta),
%       z2= sin(delta),
% for every inverter.
% No z1^2+z2^2=1 is used in the final set of equations.

xtrig = symSingleVars('xtrig',2*num_inv);
TDQ=[];
odd_aux = xtrig(1:2:end);  % odd indexes are the cosines
even_aux = xtrig(2:2:end); % even indexes are the sines

% Multilinear rotation matrix
for i=1:num_inv
    T = [odd_aux(i),-even_aux(i);even_aux(i) odd_aux(i)]; % fDQ = T*fdq
    TDQ = blkdiag(TDQ,transpose(T));
end

% Time evolution of dummy variables
eq_dz_odd = diff(odd_aux,t) == -diff(delta,t).*even_aux;
eq_dz_even = diff(even_aux,t) == diff(delta,t).*odd_aux;

% Subsitute d(delta)/dt out
eq_dz_odds = subs(eq_dz_odd,diff(delta,t),rhs(eq4));
eq_dz_evens = subs(eq_dz_even,diff(delta,t),rhs(eq4));

% These equations will have the dummy variables 
eq14s = subs(eq14,vbdq,TDQ*vbDQ);
eq15s = subs(eq15,vbdq,TDQ*vbDQ);



%% Inverter output eqns.
% Frequency
yeq1 = omega == -mp*P + wn;

% wcom is the first component of omega
eqcom = wcom == subs(rhs(yeq1(1)));

% Substitute all wcom in the equations.
eq4s = subs(eq4,wcom,rhs(eqcom));
eq10s = subs(eq10,wcom,rhs(eqcom));
eq11s = subs(eq11,wcom,rhs(eqcom));
eq12s = subs(eq12,wcom,rhs(eqcom));
eq13s = subs(eq13,wcom,rhs(eqcom));

eq14ss = subs(eq14s,wcom,rhs(eqcom)); % multilinear
eq15ss = subs(eq15s,wcom,rhs(eqcom)); % multilinear

eq14ss_nl = subs(eq14s_nl,wcom,rhs(eqcom)); % nonlinear
eq15ss_nl = subs(eq15s_nl,wcom,rhs(eqcom)); % nonlinear

% Also substitute the wcom into the dummy time evolution
eq_dz_oddss = subs(eq_dz_odds,wcom,rhs(eqcom));
eq_dz_evenss = subs(eq_dz_evens,wcom,rhs(eqcom));

% These are here because I forgot to do them above.
eq10ss = subs(eq10s,vl,rhs(eq9s));
eq11ss = subs(eq11s,vl,rhs(eq9s));

% Output DQ current variable
ioDQ = symDQvars('ioDQ',num_inv);

% Nonlinear current
yeq2_nl = ioDQ == transpose(TDQ_nl)*io;

% % Linear current
% yeq2 = ioDQ == I0DQ + Ts*(io-I0dq) + Tc*(delta-delta0);

% Multilinear current
yeq2 = ioDQ == transpose(TDQ)*io;

disp('Inverter equations in workspace')