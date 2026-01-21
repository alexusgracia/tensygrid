%% jacobian
% Calculate jacobian matrix of a dmss model.
%% Syntax
% J = jacobian(sys) Calculate jacobian matrix with default symbolic variables [xpi, xi, ui, yi].
%
% J = jacobian(sys, xp, x, u, y) Calculate jacobian matrix for given
% variable values (numeric or symbolic).
%
% J = jacobian(sys, xp, x, u, y, z) Calculate jacobian matrix for given
% variable values (numeric or symbolic) for hybrid dmss model (Experimental).
%% Description
% jacobian(sys) calculates the jacobian matrix of an implicit multilinear time-invariant (iMTI) model stored as dmss object based of the decomposed
% parameter tensor sys.H (based of its structure matrix sys.H.F and parameter matrix sys.H.phi) with numeric or smybolic parameters, using similar principles as in [1] for explicit mutlilinear time invariant (eMTI) models.
%% Examples
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
% and serves as an example system.
%
% *Calculate symbolic jacobian matrix:*
%
% The symbolic jacobian matrix can be calculated simply by
J = sys.jacobian()
%%
% For viewing the entries of the jacobian matrix the entries of the struct
% can be accessed.
disp(full([J.equality.stateDerivative, J.equality.state, J.equality.input, J.equality.algebraic]))
%%
% where the rows correspond to the six equations and the columns to the
% variables (2 state derivatives, 2 states, 4 inputs, 4 algebraic
% variables).
%
% *Calculate symbolic jacobian matrix using nondefault variable symbols:*
%
% Nondefault variable symbols can be used by using the additional inputs of
% jacobian()
dx = sym('dx',[sys.n,1]);
x = sym('x',[sys.n,1]);
in = sym('in',[sys.m,1]);
out = sym('out', [sys.p,1]);

J = sys.jacobian(dx, x, in, out);
disp(full([J.equality.stateDerivative, J.equality.state, J.equality.input, J.equality.algebraic]))

%%
% *Calculate jacobian matrix using numeric variable values:*
%
% For calculating the jacobian matrix for certain numeric variable values,
% they can be simply entered as the inputs of jacobian()
dx = zeros(sys.n, 1);
x = zeros(sys.n, 1);
u = zeros(sys.m, 1);
y = zeros(sys.p, 1);

J = sys.jacobian(dx, x, u, y);
disp(full([J.equality.stateDerivative, J.equality.state, J.equality.input, J.equality.algebraic]))

%%
% *Calculate symbolic jacobian matrix using a mix of numeric variable values and symbolic variables:*
%
% Also symbolic and numeric values can be mixed
dx = zeros(sys.n, 1);
x = sym('x',[sys.n,1]);
u = [sym("u1"); 0; 0; 0];
y = sym('y',[sys.p,1]);

J = sys.jacobian(dx, x, u, y);
disp(full([J.equality.stateDerivative, J.equality.state, J.equality.input, J.equality.algebraic]))

%%
% *Calculate jacobian matrix for numeric model parameters and variables:*
%
% Using the replaceSymbolicParameters() functions, the symbolic parameters
% of the model can be replaced by according numeric values
symParameters = [sym("m"), sym("c_d"), sym("c_v"), sym("r0"), sym("A"), sym("k")];
numParameters = [9.03*11.33*3.4*1.225, 1006, 1860, 2501*10^3, 2*(9.03+11.33)*3.4 + 9.03*11.33, 3];
sys = sys.replaceSymbolicParameters(symParameters, numParameters);

dx = zeros(sys.n, 1);
x = zeros(sys.n, 1);
u = zeros(sys.m, 1);
y = zeros(sys.p, 1);

J = sys.jacobian(dx, x, u, y);
disp(full([J.equality.stateDerivative, J.equality.state, J.equality.input, J.equality.algebraic]))

%%
% Using numeric variables as well as numeric parameters facilitates the
% fastest calculation of the jacobian.

%% References
% [1] C. Kaufmann, D. Crespí, G. Lichtenberg, Georg Pangalos, and Carlos Cateriano Yáñez, "Efficient Linearization of Explicit Multilinear Systems using Normalized Decomposed Tensors,” IFAC-PapersOnLine, vol. 56, no. 2, pp. 7312–7317, Jan. 2023, doi: https://doi.org/10.1016/j.ifacol.2023.10.344.
%% See also
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>,
% <incidenceMatrixDoc.html incidenceMatrix>
% <hyCPN1Doc.html hyCPN1>,
% <replaceSymbolicParametersDoc.html replaceSymbolicParameters>
%
%
% Author(s): Torben Warnecke
