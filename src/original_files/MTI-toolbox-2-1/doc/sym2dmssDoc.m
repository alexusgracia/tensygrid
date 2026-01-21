%% sym2dmss
% Converts symbolic or string equations into dmss object
%% Syntax
% sys = sym2dmss(eq, ts) Create a dmss object with a set of equations as
% symbolic or string array with the step size ts (ts = 0: continuous-time, ts > 0: discrete-time).
%
% sys = sym2dmss(eq, ts, symbols) Create a dmss object with a set of equations as
% symbolic or string array with the step size ts (ts = 0: continuous-time, ts > 0: discrete-time)
% with symbols as variable type symbols. The default is ["xp", "x", "u",
% "y", "z"],
% with "xp" for state predictions/derivatives, "x" for states, "u" for
% inputs, "y" for algebraic variables and "z" for binary variables
% (Experimental: used for hybrid iMTI/dmss models).
%% Description
% sym2dmss(eq,ts) creates a dmss model with an fully expanded hyCPN1
% tensor. 
%
% Eq can be a symbolic or string array (or cell) with polynomial equations 
% only containing integer exponents (auxilary variables and equations
% will be created by the algorithm to represent integer exponents not equal to 1).
% Math operators must be used in between each variable and parameter. 
% Possible operators are "*", "+", "-", "/", "^", as well as brackets "("
% and ")". For different variables of the same type use integer indices 
% after the variable type symbol, e.g. for states "x1" and "x2". 
% If an integer is missing, the leftover variables will be parameters.
% The equations do not have to be written in an implicit format ("0 =
% ..."), any format is possible, e.g. polynomials on both sides of the
% equation. Parameters can either be numeric or symbolic.
% But for simulation, symbolic values must be replaced by according numeric
% values.
%% Examples
% *Create dmss model with symbolic multilinear equations: Three-dimensional discrete-time Lotka–Volterra model*
% 
% The preditor-prey model from [1] can be described by the following
% differences equations
Eq = ["0 = xp1 - x1 * (e1 + a1_1*x1 + a1_2*x2 + a1_3*x3)";...
    "0 = xp2 - x2 * (e2 + a2_1*x1 + a2_2*x2 + a2_3*x3)";...
    "0 = xp3 - x3 * (e3 + a3_1*x1 + a3_2*x2 + a3_3*x3)"];
%%
% which can be converted into a continuous time dmss object
ts = 1;
sysD = sym2dmss(Eq, ts);
disp(sysD.symbolicEquations)
% Since default symbols have been used in the equations, no additional
% symbol array is need.
%
%%
% *Create a dmss model of a ficitve polynomial system:*
%
% For illustrative purposes a fictive polynomial system with the following
% equations will be transformed into a dmss object
Eq = ["u1*2*(dx1-x1)/(in2*out1-0.61*x1*out2^2)  = in1*in3^(-1)";...
    "out1*in1 = x1*in3";...
    "(out2-in2)^2 = x1*out1 + 0.12*in1*out1"];
%%
% Non-default symbols have been used, therefor additionally a symbol array
% is needed, indicating which symbols represent which variable types in
% order [stateDerivatives/Predictions, states, inputs, algebriac].
symbols = ["dx", "x", "in", "out"];
sys = sym2dmss(Eq, 3, symbols);

%%
% The equations have been reordered and expanded by auxiliary equations and
% variables to achieve a multilinear structure, they can be viewed by
disp(sys.symbolicEquations)

%% References
% [1] G.I. Bischi and F. Tramontana, "Three-dimensional discrete-time Lotka–Volterra models with an application to industrial clusters", Communications in Nonlinear Science and Numerical Simulation, 15(10), pp. 3000–3014., Oct. 2010. doi:10.1016/j.cnsns.2009.10.021. 

%% See also
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <replaceSymbolicParametersDoc.html replaceSymbolicParameters>,
% <hyCPN1Doc.html hyCPN1>
%
%
% Author(s): Torben Warnecke