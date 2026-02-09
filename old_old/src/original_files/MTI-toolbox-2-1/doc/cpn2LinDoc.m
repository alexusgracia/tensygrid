%% cpn2Lin 
% Computes the Jacobian of a CPN1 MTI object around an operating point
%% Syntax
%    cpn2Lin(obj,x)
%% Description
% |cpn2Lin(obj,x)| Computes the Jacobian of an explicit multilinear model stored in the CP-decomposed norm-1 format |obj| around an operating point |x|. For further information on the implementation of the algorithm see [1].     
%% Input Arguments 
% |obj|: CPN1 object 
%
% |x|  : operating point, scalar or vector 
% 
%% Output Arguments
% |J|: Jacobian matrix
%
%% Example:
% Assume we have the following thrird-order explicit MTI model $$\dot{\mathbf{x}}=\left(\begin{array}{cc} x_2 + x_1x_2\\ 2x_2 + 3x_3 + 2x_1x_2 - 3x_1x_3 \\ 4x_3 - 4x_1x_3
% \end{array}\right)$$. We write the explicit MTI model as a |CPN1| object by defining the structural matrix S to form the monomial  
%
    S=[0.5 -0.5; 1 0; 0 1];
%%
% and the parameter matrix with corresponding coefficients of the summands:
%
    phi=[2 0; 4 6; 0 8];
%%
% Then we can create explicit MTI model as a CPN1 object
%
    obj=CPN1(S,phi);
%%
% and assume the operating point
%
    x_op=[-1;1;3]; 
%%
% Then we can calcute the Jacobian
%
    J = cpn2Lin(obj,x_op)
%%
%% References
% [1] C. Kaufmann, D. Crespí, G. Lichtenberg, G. Pangalos, and C. Cateriano Yáñez, "Efficient Linearization of Explicit Multilinear Systems using Normalized Decomposed Tensors," IFAC-PapersOnLine, vol. 56, no. 2, pp. 7312–7317, Jan. 2023, doi: https://doi.org/10.1016/j.ifacol.2023.10.344.
%% See also
% <jacobianDoc.html jacobian>, <cpn2LinSparseDoc.html cpn2LinSparse>, <linearizeDoc.html linearize>