%% replaceSymbolicParameters
% Replaces symbolic parameter values of the model tensor by new (numeric or symbolic) values.
%% Syntax
% sys = sys.replaceSymbolicParameters(OLD, NEW)
% with OLD as the old symbolic parameter values as symbolic array, which will be
% replaced by new parameter values in NEW, either a symbolic or numeric array with the
% same length as OLD
%
% Attention: For simulation it is mandatory that the parameters have
% numeric values.
%% Description
% The algorithm simply replaces symbolic paramter values with new parameter values which are stored in the
% continuous sparse matrices of sys.H, 
% e.g. sys.H.F.state.c or sys.H.phi.equality.c
%% Example
% *Sprott B Attractor*
%
% In [1] J. C. Sprott proposed chaotic attractor named Sprott A-S systems. 
% The Sprott B attractor can be represented by the following model tensor
H = hyCPN1();
H.F.stateDerivative  = [1 0 0 0 0 0 0 0;...
                        0 1 0 0 0 0 0 0;...
                        0 0 1 0 0 0 0 0];
H.F.state = [0 0 0 0 1 0 0 1;...
             0 0 0 1 0 1 0 1;...
             0 0 0 1 0 0 0 0];
H.phi.equality =   [-1 0 0 sym("a") 0 0 0 0;...
                    0 -1 0 0 1 -sym("b") 0 0;...
                    0 0 -1 0 0 0 sym("c") -1];
%%
% A dmss model can be created from the tensor
Ts = 0;
sys = dmss(H,Ts);
%%
% Its equations can be visited by
disp(sys.symbolicEquations)
%%
% The symbolic parameters
params = [sym("a"), sym("b"), sym("c")];
%% 
% can be replaced by numeric values
numParams = [0.4, 1.2, 1];
%%
% using
sys = sys.replaceSymbolicParameters(params, numParams);
disp(sys.symbolicEquations)
%% References
% [1] J. C. Sprott, “Some simple chaotic flows,” Physical Review E, vol. 50, no. 2, Aug. 1994. doi:10.1103/physreve.50.r647 
%% See also
% <dmssDoc.html dmss>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>
%
%
% Author(s): Torben Warnecke