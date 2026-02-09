%% d2c
% Converts discrete-time dynamic iMTI system to continuous-time.

%% Syntax
% SYSC = d2c(SYSD,TS) computes a continuous-time model with using one step
% forward euler approximation.
%
%% Description
% For the continuization the one step forward euler approximation is used: xd = (xp - x)/TS.
% The algorithm directly implements additional auxilary equations for all state derivatives 
% and replaces the state predicitons by auxiliary variables by manipulation
% of the parameter tensor (the strucutre and parameter matrix of the CP
% decomposition).
% Afterwards it will be checked if those additional auxilary variables and
% equations can be reduced algebraicly, if so, they will be eliminated by
% algbraic elimination by direct manipulation of the parameter tensor 
% (the strucutre and parameter matrix of the CP decomposition).
%
%% Examples
% *Calculate continuous-time model of a 3 dimensional discrete-time lotka-volterra model:*
%
% The preditor-prey model from [1] can be described by the following
% differences equations
Eq = ["xp1 = x1 * (e1 + a1_1*x1 + a1_2*x2 + a1_3*x3)";...
    "xp2 = x2 * (e2 + a2_1*x1 + a2_2*x2 + a2_3*x3)";...
    "xp3 = x3 * (e3 + a3_1*x1 + a3_2*x2 + a3_3*x3)"];
%%
% which can be converted into a continuous time dmss object
ts = 1;
sysD = sym2dmss(Eq, ts);
disp(sysD.symbolicEquations)
%%
% For simulation the symbolic parameters should be replace by numeric
% values
e = sym(["e1", "e2", "e3"]);
a = reshape(sym("a", [3,3]), 1,[]);
symParameters = [e, a];
numParameters = [2, 2, 2, -1, -0.5, 0.61, 0.61, -1, -0.5, -0.5, 0.61, -1];
sysD = sysD.replaceSymbolicParameters(symParameters, numParameters);
disp(sysD.symbolicEquations)

%%
% The discrete-time system can be converted into a continuous-time system
% by
sysC = sysD.d2c();
disp(sysC.symbolicEquations)

%%
% Simulating the two system, their behaviour can be compared
x0 = [0.8 0.8 1.5];
tend = 1e2;
simoutD = dmsim(sysD, x0, [0:ts:tend]);

simoutC = dmsim(sysC, x0, [0, tend]);

figure()
plot3(simoutC.x(:,1),simoutC.x(:,2),simoutC.x(:,3), "r-", LineWidth=0.01);
hold on
plot3(simoutD.x(:,1),simoutD.x(:,2),simoutD.x(:,3), "k--", LineWidth=0.01);
view([130, 15])
grid on
legend(["continuous model", "discrete model"], "Location","southeast")
xlabel('x1')
ylabel('x2')
zlabel('x3')
hold off

figure()
plot(simoutC.tsim, nan(length(simoutC.tsim),1), "r-", LineWidth=0.01);
hold on
plot(simoutD.tsim, nan(length(simoutD.tsim),1), "k--", LineWidth=0.01);
plot(simoutC.tsim, simoutC.x, "r-", LineWidth=0.01);
plot(simoutD.tsim, simoutD.x, "k--", LineWidth=0.01);
grid on
legend(["continuous model", "discrete model"], "Location","southeast")
xlabel('time')
ylabel('x')
hold off


%%
% While the continuous-time system reaches a stable non-periodic
% trajectory, the discrete-time model reaches a 12-cycle.

%% References
% [1] G.I. Bischi and F. Tramontana, "Three-dimensional discrete-time Lotka–Volterra models with an application to industrial clusters", Communications in Nonlinear Science and Numerical Simulation, 15(10), pp. 3000–3014., Oct. 2010. doi:10.1016/j.cnsns.2009.10.021. 

%% See also
% <dmssc2dDoc.html c2d>,
% <d2dDoc.html d2d>,
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>
%
%
% Author(s): Torben Warnecke