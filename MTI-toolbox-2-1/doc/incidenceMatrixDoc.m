%% incidenceMatrix
% Calculate the incidence matrix of a dmss model.
%% Syntax
% I = incidenceMatrix(sys) Calculates the incidence matrix I of a dmss model
% sys, directly from the parameter tensor sys.H. I is stored as struct with entries for equations (and inequalities for hybrid dmss models (experimental)).
% Additionally it stores entries for each variable type.
%% Example
% The following equations describe a simple HVAC system
Eq =                            ["m*xp1 - y4 - y2 = 0";...
                                      "m*xp2 - y3 = 0";...
                 "x1 - c_d*y1 - r0*x2 - c_v*x2*y1 = 0";...
"y2 + u1*x1 - c_d*u1*u2 - r0*u1*u3 - c_v*u1*u2*u3 = 0";...
                              "y3 - u1*u3 + u1*x2 = 0";...
                            "y4 + A*k*u4 - A*k*y1 = 0"];
%%
% which can be converted into a continuous-time dmss model by
sys = sym2dmss(Eq, 0);

%%
% The incidence matrix can be calculated by
I = incidenceMatrix(sys)

%%
% For viewing the entries of the incidence matrix the entries of the struct
% can be accessed.
disp(full([I.equality.stateDerivative, I.equality.state, I.equality.input, I.equality.algebraic]))
%%
% where the rows correspond to the six equations and the columns to the
% variables (2 state derivatives, 2 states, 4 inputs, 4 algebraic
% variables). Entries of 1 and 0 represent whether a variable is occuring
% in the according equation or not.

%% References
% [1] ...CoDIT2024...
%% See also
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>,
% <dmssJacobianDoc.html jacobian>
%
%
% Author(s): Torben Warnecke