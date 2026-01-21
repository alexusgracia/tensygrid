% This scripts puts the inverter and network-load equations together.

%% Connecting inverters and network.
% Network Voltage nodes are equal to the bus voltages of the inverters in DQ.
eq14sss = subs(eq14ss,vbDQ,Agen'*v);
eq15sss = subs(eq15ss,vbDQ,Agen'*v);

eq14sss_nl = subs(eq14ss_nl,vbDQ,Agen'*v);
eq15sss_nl = subs(eq15ss_nl,vbDQ,Agen'*v);

% Generated currents in the network are the ouput current of the inverters in DQ.
eq2bls = subs(eq2bl,Igen,rhs(yeq2));
eq2bls_nl = subs(eq2bl,Igen,rhs(yeq2_nl));

% Subsitute wcom in the network equations.
eq1bls = subs(eq1bl,wcom,rhs(eqcom));

%% Multilinear Equations with linearized trigs.
% disp("Linearized Trigs for multilinear equations ");
% % Collection of all equations in the inverters and network.
% % The final "version" of an equations is the one with the more "s".
% eqns = [eq1;eq2;eq3;eq4s;eq5;eq6;eq7s;eq8s;eq10ss;eq11ss;eq12s;eq13s;eq14sss;eq15sss;eq1bls;eq2bls];
% eqns = eqns(t);
% % State vector, these are the only variables that should be in the equations
% % that's why there were so many substitutions.
% vars = [idc;vdc;delta;P;Q;phi;gamma;il;vo;io;f;v;zdc];
% disp("Number of Variables: "); origVars = length(vars)
% disp("Number of Equations: "); origeqns = length(eqns)
% disp("Eqns check")

%% Multilinear Equations with dummy xtrig and no sin^2+cos^2 = 1 constraints.

% Collect all equations
eqns = [eq1;eq2;eq3;eq4s;eq5;eq6;eq7s;eq8s;eq10ss;eq11ss;eq12s;eq13s;eq14sss;eq15sss;eq1bls;eq2bls;eq_dz_oddss;eq_dz_evenss];
eqns = eqns(t);
% Variable vector, these are the only variables that should be in the equations
% vars = [idc;vdc;delta;P;Q;phi;gamma;il;vo;io;f;v;zdc;xtrig]; 

% Order is vars = [diff. variables, algebraic variables] 
vars = [idc;vdc;delta;P;Q;phi;gamma;il;vo;io;f;xtrig;v;zdc]; 



%% Nonlinear Equations 
eqns_nl = [eq1;eq2_nl;eq4s;eq5;eq6;eq7s;eq8s;eq10ss;eq11ss;eq12s;eq13s;eq14sss_nl;eq15sss_nl;eq1bls;eq2bls_nl];
eqns_nl = eqns_nl(t);

vars_nl = [idc;vdc;delta;P;Q;phi;gamma;il;vo;io;f;v];



