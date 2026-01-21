%% hyCPN1
% hybrid CP decomposed norm-1 tensor for dmss model parameters
%% Syntax
% H = hyCPN1() Create empty hyCPN1 tensor with hybrid parameter structure.
% For information on how to access and change model parameters look into examples.
%
%% Description
% hyCPN() stores the model paramters of a dmss model in the structure of a hybrid norm-1 CP
% decomposed tensor.
% The norm-1 CP-decomposed tensor format has been introduced in [2].
% A hybrid implicit multilinear time invariant (hybrid iMTI) model has been
% briefly introduced in [1].
%
%
% *hybrid structure:* 
%
% For the structure matrix H.F, the rows of H.F have been split into the
% different variable types:
%
% H.F.state = [];
%
% H.F.stateDerivative = [];
% 
% H.F.algebraic = [];
% 
% H.F.boolean = []; --> Experimental: For hybrid dmss models
% 
% H.F.input = [];
% 
%
% For the parameter matrix H.phi, the rows of H.phi have been split into the
% different equations and inequality constraints (for hybrid dmss models):
%
% H.phi.equality = [];
% 
% H.phi.inequality = []; --> Experimental: For hybrid dmss models
%
%
% *Matrix structure:*
% 
% Matrices in these structures will stored in a hybrid fashion with the
% following sub structure as 
%
% For the structure matrix H.F:
%
% -> H.F.state.t stores all true*x (equivalent to 1*x) values of  H.F.state
% as sparse logical matrix
%
% -> H.F.state.f stores all false*x (equivalent to (1-x)) values of
% H.F.state as sparse logical matrix
%
% -> H.F.state.c stores all continuous valued parameters (equivalent to
% (1-F + F*x)) of H.F.state as sparse double/symbolic matrix
%
%
% For the parameter matrix H.phi:
%
% -> H.phi.equality.t stores all true (equivalent to 1) values of
% H.F.state as sparse logical matrix
%
% -> H.F.state.f stores all false (equivalent to -1) values of  H.F.state as sparse logical matrix
%
% -> H.F.state.c stores all continuous valued parameters of H.F.state as sparse double/symbolic matrix
%
%% Examples
% To generate a model using an hyCPN1 object, first an empty hyCPN1 object
% must be created
H = hyCPN1()

%%
% Now the structure of H can be filled with matrices by accessing its
% properties. It can be either filled by logical, double or symbolic
% values.
H.F.stateDerivative = [1 0 0 0 0 0 0 0];
H.F.state = [0 1 0 0 1 0 sym("k_3") sym("k_5")];
H.F.input = [0 1 1 0 0 0 sym("k_2") sym("k_4");...
            0 0 1 0 0 0 0 0;...
            0 0 0 0 0 1 0 0];
H.F.algebraic = [0 0 0 1 0 0 0 0;...
                0 0 0 0 1 1 sym("k_1") 0];

H.phi.equality = [sym("m_w")*sym("c_w") sym("c_w") -sym("c_w") 1 0 0 0 0;...
                    0 0 0 1 -sym("A") sym("A") 0 0;...
                    0 0 0 0 0 0 sym("phi_1") sym("phi_2")];
%%
% When setting a sub matrix (e.g. H.F.state) the matrix will be stored in a
% sparse hybrid fashion (H.F.state.t, H.F.state.f and H.F.state.c), for making use
% of faster sub-routines in simulation etc. Moreover these metrices are
% used in different routines, like in connecting different dmss models
% witch connect().
% Also those sub matrices can be
% accessed, set or changed directly, e.g. via H.F.state.t = [0 1 0 0 1 0 0
% 0].
%
% Now a dmss object can be created from the filled hyCPN1 object.
sys = dmss(H)

%%
% Symbolic equations can be generated with
disp(sys.symbolicEquations)
%%
% e.g. for checking the correctness model. Symbolic parameters can be
% changed by
symParameters = [sym("m_w") sym("c_w") sym("A") sym("phi_1") sym("phi_2") sym("k_1") sym("k_2") sym("k_3") sym("k_4") sym("k_5")];

di = 206.5/1000;
A = di^2 *pi()/4;
m_w = A*998*1; 
c_w = 4200;
k = [-0.1497, -2.3659, -0.0109, 0.5556, 0.0131, 0.9872, 0.0259];
numParameters = [m_w, c_w, A, k];

sys = sys.replaceSymbolicParameters(symParameters, numParameters);

symParameters = [sym("m_w") sym("c_w") sym("A") sym("phi_1") sym("phi_2") sym("k_1") sym("k_2") sym("k_3") sym("k_4") sym("k_5")];

disp(sys.symbolicEquations)

%% References
% [1] G. Lichtenberg et al., "Implicit multilinear modeling," at - Automatisierungstechnik, vol. 70, no. 1, pp. 13–30, Jan. 2022. doi:10.1515/auto-2021-0133 
%
% [2] N. Jöres et al., “Reduced CP representation of multilinear models,” Proceedings of the 12th International Conference on Simulation and Modeling Methodologies, Technologies and Applications, 2022. doi:10.5220/0011273100003274  
%
%% See Also
% <dmssDoc.html dmss>,
% <sym2dmssDoc.html sym2dmss>,
% <mss2dmssDoc.html mss2dmss>,
% <connectDoc.html connect>,
% <replaceSymbolicParametersDoc.html replaceSymbolicParameters>
%
%
% Author(s): Torben Warnecke