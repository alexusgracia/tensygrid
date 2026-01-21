%% dmss
% Multilinear descriptor state-space model (implicit)
%% Syntax
% sys = dmss() 
% Create empty dmss object with an empty hyCPN1 tensor/object.
%
% sys = dmss(H) 
% Construct a continuous-time (ts = 0) dmss object 
% when H is a hyCPN1 tensor/object.
%
% sys = dmss(mssSys) 
% Convert a mss object into a dmss object with a similar 
% stepsize ts, when mssSys is a mss object.
%
% sys = dmss(H,ts) 
% Construct a dmss object when H is a hyCPN1 objct/tensor, 
% with ts as stepsize (ts = 0: continuous-time, ts > 0: discrete-time).

%% Description
% The implicit multilinear time invariant (iMTI) model format has been introduced
% in [1]. In general, iMTI models can be described by a set of $N$ implicit continuous-
% time differential algebraic equations (DAEs) 
%
% $$ 0 = \mathbf{h}(\mathbf{\dot{x}}, \mathbf{x}, \mathbf{u}, \mathbf{y}) $$
%
% containing only multilinear functions $\mathbf{h}$.
% Multilinear models can be represented using tensors and the system of equations can be represented by using the extended tensor
% product 
%
% $$0 = \langle \mathbf{H} \vert \mathbf{M}(\mathbf{\dot{x}}, \mathbf{x}, \mathbf{u}, \mathbf{y}) \rangle$$ 
%
% of the model tensor $\mathbf{H}$ with the dimensions $n\times n\times m\times p\times N$ and monomial tensor $\mathbf{M}$ with the dimensions $n\times n\times m\times p$, see [1].
% This utilizes tensor decomposition methods which introduces the reduced norm-1
% canonical polyadic decomposed (CPn1) tensor representation of iMTI
% models, see [2],
% where the $j^{th}$ equation 
%
% $$ 0 = \sum_{k=1}^{R} \Phi_{j,k} \prod_{l=1}^{2n+m+p} \left(1 - \vert F_{l,k} \vert + F_{l,k} v_{l}\right)$$
% 
% can be represented in terms of the structure matrix $F$ with the dimensions $(2n+m+p)\times R$ and parameter
% Matrix $\Phi$ with the dimensions $N\times R$. While
% 
% $$\mathbf{v} = \left[\mathbf{\dot{x}}, \mathbf{x}, \mathbf{u}, \mathbf{y}\right]'$$ 
% 
% represents the signal vector of state derivatives (continuous time) / next states (discrete time), inputs and auxiliary/algebriac variables.
%
% sys = dmss(H, ts) creates an implicit/decriptor multilinear time invriant 
% (iMTI) model stored by using the sparse norm-1 canonical polyadic 
% decomposed (CPn1) tensor H, as a hyCPN1 object, 
% and the stepsize ts.
% 
%% Properties
% <html><h3>H - System tensor</h3></html> 
%
% System tensor, specified as canonical decomposed norm-1 tensor
% (hyCPN1-object).
%
% <html><h3>n - Number of states</h3></html> 
%
% Number of states, specified as an integer $n\geq 0$.
%
% <html><h3>m - Number of inputs</h3></html> 
%
% Number of inputs, specified as an integer $m\geq 0$.
%
% <html><h3>p - Number of algebraic variables</h3></html> 
%
% Number of states, specified as an integer $q\geq 0$.
%
% <html><h3>nEq - Number of equations</h3></html>
%
% Number of equations, specified as an integer $nEq\geq 0$.
%
% <html><h3>stateName - Names of states</h3></html> 
%
% Names of states, specified as an string vector with $n$ entries.
%
% <html><h3>stateUnit - Units of states</h3></html>
%
% Units of states, specified as an string vector with $n$ entries.
%
% <html><h3>algebraicName - Names of algebraic variables</h3></html>
%
% Names of algebraic variables, specified as an string vector with $p$ entries.
%
% <html><h3>algebraicUnit - Units of algebraic variables</h3></html>
%
% Units of algebraic variables, specified as an string vector with $p$ entries.
%
% <html><h3>inputName - Names of inputs</h3></html> 
%
% Names of inputs, specified as an string vector with $m$ entries.
%
% <html><h3>inputUnit - Units of inputs</h3></html>
%
% Units of inputs, specified as an string vector with $m$ entries.
%
%% Object Functions
%
% <html>
% <table border=1>
% <tr>
% <td>dmsim</td>
% <td>Simulation class to simulate time response of dynamic system to arbitrary inputs.</td>
% </tr>
% <tr>
% <td>sym2dmss</td>
% <td>Transfer symbolic equations into dmss model.</td>
% </tr>
% <tr>
% <td>c2d</td>
% <td>Convert continuous-time dynamic iMTI system to discrete time.</td>
% </tr>
% <tr>
% <td>d2c</td>
% <td>Convert discrete-time dynamic iMTI system to continuous-time.</td>
% </tr>
% <tr>
% <td>d2d</td>
% <td>Resample discrete-time dynamic iMTI system.</td>
% </tr>
% <tr>
% <td>append</td>
% <td>Append multiple dmss-objects.</td>
% </tr>
% <tr>
% <td>connect</td>
% <td>Connect signals of dmss objects, based of index matrices or signal names.</td>
% </tr>
% <tr>
% <td>jacobian</td>
% <td>Calculate jacobian matrix of the set of equations of a dmss model.</td>
% </tr>
% <tr>
% <td>incidenceMatrix</td>
% <td>Calculate incidence matrix of the set of equations of a dmss model.</td>
% </tr>
% <tr>
% <td>symbolicEquations</td>
% <td>Calculate symbolic equations of a dmss model with symbolic variables</td>
% </tr>
% <tr>
% <td>algebraicElimination</td>
% <td>Eliminate algrbaic equalities and variables of a dmss model, if its possible while maintaining the implicit multilinear structure.</td>
% </tr>
% <tr>
% <td>normalize</td>
% <td>Coordinate transformation of variables for dmss models.</td>
% </tr>
% <tr>
% <td>replaceSymbolicParameters</td>
% <td>Replace symbolic parameter values of the model tensor.</td>
% </tr>
% </table>
% </html>
%% Examples
% *Symbolic to dmss conversion and simulation of the the aizawa attractor*
%
% The aizawa attractor (from [3]) can be represented by an iMTI model using
% additional auxiliary variables and euqations. When using sym2dmss needed
% auxiliary variables will be generated automatically.
% The equations of the aizawa attractor are
Eq = {"dx1 = (x3-b)*x1 -d*x2";...
    "dx2 = d*x1+ (x3-b)*x2";...
    "dx3 = c+ a*x3 - x3^3/3 - x1^2 + f*x3*x1^3"};
%%
% and can be directly used for creating a continuous time dmss object with
Ts = 0;
sys = sym2dmss(Eq, Ts, {"dx", "x", "u", "y"});
%%
% The equations of the resulting iMTI model can be viewed with
symEq = sys.symbolicEquations
%%
% where $xp_i$ are the state derivatives, $x_i$ are states (differential
% variables), $y_i$ are algebriac variables and $u_i$ are input signals.
% One can see that due to the quadratic and cubic terms additional
% auxiliary variables y along with additional auxiliary equations have been created.
%
% For simulation the symbolic Parameters must be replaced by numeric values
symParameters = ["a", "b", "c", "d", "f"];
numParameters = [0.095, 0.7, 0.65, 3.5, 0.1];
sys = sys.replaceSymbolicParameters(symParameters, numParameters);
%%
% Simulation can be done with creating a dmsim object (for more see: ...) using
% the iMTI model sys
tend = 1e3;
x0 = [0.1, 0, 0];
simout = dmsim(sys, x0, [0 tend], []);
%%
% The results can be accessed as properties of the dmsim object and plotted
plot3(simout.x(:,1), simout.x(:,2), simout.x(:,3))
grid on

%%
% *Creation by using the hyCPN1 tensor class for a simple thermal flow system*
%
% By creating a empty hyCPN1 object the object structure will be generated
H = hyCPN1();
%%
% and can be accessed via the properties of the object.
% The structure matrix can be manipulated by calling its property
H.F.stateDerivative = [1 0 0 0 0 0 0 0];
H.F.state = [0 1 0 0 1 0 sym("k_3") sym("k_5")];
H.F.input = [0 1 1 0 0 0 sym("k_2") sym("k_4");...
            0 0 1 0 0 0 0 0;...
            0 0 0 0 0 1 0 0];
H.F.algebraic = [0 0 0 1 0 0 0 0;...
                0 0 0 0 1 1 sym("k_1") 0];
%%
% Similar the parameter matrix can be accessed, here using partly symbolic
% parameters.
H.phi.equality = [sym("m_w")*sym("c_w") sym("c_w") -sym("c_w") 1 0 0 0 0;...
                    0 0 0 1 -sym("A") sym("A") 0 0;...
                    0 0 0 0 0 0 sym("phi_1") sym("phi_2")];
%%
% A continuous time dmss model can be created by
sys = dmss(H,0);
%%
% The equations can be checked by
disp(sys.symbolicEquations)
%%
% where $xp_i$ are the state derivatives, $x_i$ are states (differential
% variables), $y_i$ are algebriac variables and $u_i$ are input signals.
%
% As before the symbolic parameters can be changed to numeric parameters by
% the replaceSymbolicParmeters function.
symParameters = [sym("m_w") sym("c_w") sym("A") sym("phi_1") sym("phi_2") sym("k_1") sym("k_2") sym("k_3") sym("k_4") sym("k_5")];

di = 206.5/1000;
A = di^2 *pi()/4;
m_w = A*998*1; 
c_w = 4200;
k = [-0.1497, -2.3659, -0.0109, 0.5556, 0.0131, 0.9872, 0.0259];
numParameters = [m_w, c_w, A, k];

sys = sys.replaceSymbolicParameters(symParameters, numParameters);

disp(sys.symbolicEquations)
%%
% The system can be simulated by using the dmsim class.
% A time vector, initial states and input values need to be defined.
t = 0:60;
x0 = 8;
u = ones([length(t),3]);
u(:,1) = 0.8*A*998 * (0.6 + 0.2*rand(length(t),1));
u(:,2) = 8 + rand(1,length(t));
u(:,3) = 25 + asinh((t-30)/10).*cos(t/5);
%%
% For continuous time system standard ode options can be used
opt = odeset('RelTol',1e-5,'AbsTol', 1e-3, 'MaxStep', 0.5);
%%
% and simulation will be performed
simout = dmsim(sys, x0, t, u, opt);

%%
% The results can be accessed and plotted by the properties of the
% simulation object simout
subplot(2,1,1)
plot(simout.tsim, simout.x)
grid on
subplot(2,1,2)
plot(simout.ut, simout.u)
grid on
%% References
% [1] G. Lichtenberg et al., "Implicit multilinear modeling," at - Automatisierungstechnik, vol. 70, no. 1, pp. 13–30, Jan. 2022. doi:10.1515/auto-2021-0133 
%
% [2] N. Jöres et al., “Reduced CP representation of multilinear models,” Proceedings of the 12th International Conference on Simulation and Modeling Methodologies, Technologies and Applications, 2022. doi:10.5220/0011273100003274  
%
% [3] W. F. Langford, "Numerical studies of torus bifurcations,” Numerical Methods for Bifurcation Problems, pp. 285–295, 1984. doi:10.1007/978-3-0348-6256-1_19  
%% See Also
% <mssDoc.html mss>,
% <msimDoc.html msim>,
% <dmsimDoc.html dmsim>,
% <sym2dmssDoc.html sym2dmss>,
% <mss2dmssDoc.html mss2dmss>,
% <dmssJacobianDoc.html jacobian>,
% <incidenceMatrixDoc.html incidenceMatrix>,
% <symbolicEquationsDoc.html symbolicEquations>,
% <dmssc2dDoc.html c2d>,
% <d2cDoc.html d2c>,
% <d2dDoc.html d2d>,
% <algebraicEliminationDoc.html algebraicElimination>,
% <normalizeDoc.html normalize>,
% <replaceSymbolicParametersDoc.html replaceSymbolicParameters>,
% <appendDoc.html append>,
% <connectDoc.html connect>,
% <hyCPN1Doc.html hyCPN1>
%
%
% Author(s): Torben Warnecke
