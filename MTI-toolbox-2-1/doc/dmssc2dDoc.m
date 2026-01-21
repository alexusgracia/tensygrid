%% c2d
% Convert model from continuous to discrete time
%% Syntax
% SYSD = c2d(SYSC,TS) computes a discrete-time model SYSD with time step TS using one step
%   forward euler approximation from an continuous-time model SYSC.
%% Description
% For the discretization the one step forward euler approximation is used: xd = (xp - x)/TS.
% The algorithm directly implements additional auxilary equations for all state predictions 
% and replaces the state derivatives by auxiliary variables by manipulation
% of the parameter tensor (the strucutre and parameter matrix of the CP
% decomposition).
% Afterwards it will be checked if those additional auxilary variables and
% equations can be reduced algebraicly, if so, they will be eliminated by
% algbraic elimination by direct manipulation of the parameter tensor 
% (the strucutre and parameter matrix of the CP decomposition).
%% Examples
% *Discretization and simulation of the Dadras Attractor using the
% dmss-class and c2d function*
%
% The Dadras attractor, introduced in [1], can be described by the
% following set of ODEs
eq = ["xp1 = x2 - a*x1 + b*x2*x3";...
    "xp2 = c*x2 - x1*x3 + x3";...
    "xp3 = d*x1*x2 - h*x3"];
%%
% which can be converted into a continuous time dmss object
ts = 0;
sys = sym2dmss(eq, ts);
disp(sys.symbolicEquations)
%%
% For simulation the symbolic parameters should be replace by numeric
% values
symParameters = sym(["a", "b", "c", "d", "h"]);
numParameters = [3, 2.7, 4.7, 2, 9];
sys = sys.replaceSymbolicParameters(symParameters, numParameters);
disp(sys.symbolicEquations)
%%
% The model can be converted into discrete time model by
tsD = 0.01;
sysD = c2d(sys, tsD);
disp(sysD.symbolicEquations)

%%
% With simulation the phase potrait can be computed for the continuous- and discrete-time models. Also the simulation times can be
% compared.

x0 = [5, 0, -4];
tend = 1e2;

tic()
simoutC = dmsim(sys, x0, [0, tend]);
timeC = toc();

tic()
simoutD = dmsim(sysD, x0, [0:tsD:tend]);
timeD = toc();

figure()
plot3(simoutC.x(:,1),simoutC.x(:,2),simoutC.x(:,3), "r-", LineWidth=0.01);
hold on
plot3(simoutD.x(:,1),simoutD.x(:,2),simoutD.x(:,3), "k--", LineWidth=0.01);
view([45, 30])
grid on
legend([sprintf("continuous model: %.3f sec.", timeC), sprintf("discrete model: %.3f sec.", timeD)], "Location","southeast")
hold off
%%
% In the phase potrait it can be seen that the discretization is stable and
% behaves similar (within a similar magnitude) as the continuous time
% model. Due to the chaotic behaviour of the system it is extremely unlikely to
% generate exactly similar trends.

%% References
% [1] S. Dadras and H. R. Momeni, “A novel three-dimensional autonomous chaotic system generating two, three and four-scroll attractors,” Physics Letters A, vol. 373, no. 40, pp. 3637–3642, Sep. 2009. doi:10.1016/j.physleta.2009.07.088
%
%% See also
% <d2cDoc.html d2c>,
% <d2dDoc.html d2d>,
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>
%
%
% Author(s): Torben Warnecke
