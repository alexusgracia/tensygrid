%% linearize
% Linearize an explicit (mss) or hybrid implicit (dmss) multilinear model around an operating point.  
%% Syntax
%   lsys = linearize(sys, op)                  % mss -> ss/sparss
%   dss  = linearize(sys, op)                  % dmss -> dss
%   [dss, dssSparse] = linearize(sys, op)      % dmss -> dense and sparse
%
% Backward compatible signatures:
%   linearize(msys, x, u)                      % mss (legacy)
%   linearize(msys, xp, x, u, y, z)            % dmss (legacy)
%% Description
% |linearize(sys, op)| invokes the class-specific linearization implemented by
% both |mss| and |dmss| objects. It accepts a single operating point argument |op|
% as a struct for clarity and extensibility.
%
% For |mss| objects:
%
% * |op.x| (n×1) – state operating point (required)
% * |op.u| (m×1) – input operating point (required if m>0)
%
% For |dmss| objects:
%
% * |op.xp| (n×1) – state-derivative operating point (optional)
% * |op.x|  (n×1) – state operating point (optional)
% * |op.u|  (m×1) – input operating point (optional)
% * |op.y|  (p×1) – algebraic variable operating point (optional)
% * |op.z|  (q×1) – boolean variable operating point (optional)
%
% Any omitted fields default to |[]|. The class methods perform basic dimension checks.
%
% Convenience forms are supported for parity with MATLAB’s style:
%
% * |mss|: numeric vector |[x; u]| or cell |{x, u}|
% * |dmss|: numeric vector |[xp; x; u; y; z]| or cell |{xp, x, u, y, z}|
%% Input Arguments
% |sys|: multilinear MTI model (|mss| or |dmss| object)
%
% |op| : struct, numeric vector, or cell containing operating point values (see above)
%% Output Arguments 
% |lsys|: linear LTI model, |ss| or |sparss| (for |mss| inputs)
%
% |dss|, |dssSparse|: linear descriptor model(s) as |dss| and |sparss| (for |dmss| inputs)

%% Example (mss):
% Second-order explicit MTI model with 2 states and 1 input (expanded and factored form)
% 
% $$\dot{\mathbf{x}}=\left(\begin{array}{cc} x_2 + x_1x_2\\ 3u_1 + 2x_2 - 3u_1x_1 + 2x_1x_2
% \end{array}\right) = \left(\begin{array}{cc} 2\cdot(0.5+0.5x_1)x_2 \\
% 4\cdot(0.5+0.5x_1)x_2 + 6 \cdot(0.5-0.5x_1)u_1 \end{array}\right) $$. 
%
% The factored explicit MTI model as |CPN1| object thus has a structure matrix  
%
    S = [0.5 -0.5; 1 0; 0 1];
%%
% and the parameter matrix
%
    phi = [2 0; 4 6];
%%
% Then we can create explicit MTI model as a |CPN1| object and the
% |mss|-object
%
    tens = CPN1(S,phi);
    msys = mss(tens);

%%
% and assume the operating point for the states |x| and the input |u|
%
    x_op = [-1;1]; 
    u_op = 3;
%%
% Then we can linearize the explicit MTI model using an op struct
%
    op = struct('x', x_op, 'u', u_op);
    lsys = linearize(msys, op)

%% Example (dmss):
% For the following continuous-time implicit multilinear model with 2 states and 1
% algebraic variable, a |dmss| object is created using the
% |sym2dmss| command:
%
    eqs = {'0==xp1-x1*x2*y1+4';
           '0==xp2-x2+2'
           '0==x2-y1';};
    dsys = sym2dmss(eqs, 0);
%%
% The operating point point, where xp=0, is: 
%
    xpop = [0; 0];
    xop  = [1; 2];
    yop  = 2;
%%
% which is stored as a |op| struct
%
    op = struct('xp', xpop, 'x', xop,'y',yop);
%%
% Then the |dmss| is linearized at the operating point as follows:
%
    dss = linearize(dsys, op);
%%
% The result |dss| is a linear descriptor state-space model, which is a
% |ss| object. 

%% References
% [1] C. Kaufmann, D. Crespí, G. Lichtenberg, G. Pangalos, and C. Cateriano Yáñez, "Efficient Linearization of Explicit Multilinear Systems using Normalized Decomposed Tensors," IFAC-PapersOnLine, vol. 56, no. 2, pp. 7312–7317, Jan. 2023, doi: https://doi.org/10.1016/j.ifacol.2023.10.344.
%
% [2] C. Kaufmann, G. Pangalos, G. Lichtenberg, O. Gomis-Bellmunt,  "Small-Signal Stability Analysis of Power Systems by Implicit Multilinear Models," preprint available on arXiv, 2025. 

%%
% You can find a full example here: <matlab:open("Demo_Linearization.mlx") open linearization example>
%% See also 
% <jacobianDoc.html jacobian>, <cpn2LinDoc.html cpn2Lin>, <cpn2LinSparseDoc.html cpn2LinSparse>