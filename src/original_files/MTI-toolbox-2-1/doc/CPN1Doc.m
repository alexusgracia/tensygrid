%% CPN1
% Canonical polyadic norm-1 form

%% Description
% Use |CPN1| to create a canonical polyadic (CP) norm-1 tensor, i.e., a CPN1 tensor.
% 
% A tensor $\mathrm{T}\in\mathbf{R}^{I_1\times I_2\times \ldots \times I_z}$ can be decomposed into a CP tensor by the sum of $r$ outer products, where $r$ is
% the so called _tensor rank_. All factors are represented in the so called factor 
% matrices $\mathbf{F}_i\in\mathbf{R}^{I_i\times r}$ for each
% dimension $i$ of size $I_i$ of the original tensor $\mathrm{T}$ up to $z$, which can be
% abreviated as
%
% $$
% \mathrm{T} = \left[\mathbf{F}_1, \mathbf{F}_2, \ldots,
% \mathbf{F}_z\right].
% $$
%
% An element of $\mathrm{T}$ is then given by
%
% $$
% \mathrm{T}\left(i_1, i_2, \ldots, i_z\right) = \sum_{j=1}^{r}
% \mathbf{F}_1\left(i_1, j\right)\cdot\mathbf{F}_2\left(i_2,
% j\right)\cdot\ldots\cdot\mathbf{F}_z\left(i_z, j\right) .
% $$
%
% A CPN1 tensor is a CP tensor where
%
% $$
% \Vert\mathbf{F}_i\left( I_i, k\right)\Vert_1 = 1 ,
% $$
% for all $k$ from $1$ to $r$ and for each dimension $i$.
%
% Particularly for mti systems, the parameter tensor factor matrices have dimensions 
% $\mathbf{F}_i\in\mathbf{R}^{2\times r}$ for the first $z-1$ dimensions
% corresponding to the monomial dimensions and the last $z$ dimension called $\phi$ is saved for
% parameterization, i.e., $\mathbf{F}_{\phi}\in\mathbf{R}^{I_z\times r}$,
% such that
%
% $$
% \mathrm{T} = \left[\mathbf{F}_1, \mathbf{F}_2, \ldots,
% \mathbf{F}_{z-1}, \mathbf{F}_{\phi}\right].
% $$
%
% For further reference see [1].
%
%% Creation
% <html><h3>Syntax</h3></html> 
%
% |CPN1tens = CPN1(Umat,Phimat)|
%
% |CPN1tens = CPN1(tens)| 
%
% <html><h2>Description</h2></html>
%
% |CPN1tens = CPN1(Umat,Phimat)| creates a CP norm-1 tensor |Ttens| of the
% form
% $$
% \mathrm{T} = \left[\mathbf{F}_1, \mathbf{F}_2, \ldots,
% \mathbf{F}_{z-1}, \mathbf{F}_{\phi}\right].
% $$
%
% |CPN1tens = CPN1(tens)| takes a CP decomposed _ktensor_ of the form
% $$
% \mathrm{T} = \left[\mathbf{F}_1, \mathbf{F}_2, \ldots,
% \mathbf{F}_{z-1}, \mathbf{F}_{\phi}\right]
% $$
%  from the Tensor Toolbox and creates a CPN1 tensor object 
%
% <html><h2>Input Arguments</h2></html>
%
% <html><h3>Umat - Structure matrix</h3></html> 
%
% Structure matrix, which in the context of the MTI, vertically stacks the
% second row of each CPN1 factor matrix corresponding to each dimension of
% the monomial. For example for explicit MTI systems:
%
% $$
% \mathbf{U}=\left[\begin{array}{c}\mathbf{F}_{x_1}(2,:)\\\mathbf{F}_{x_2}(2,:)\\\vdots\\\mathbf{F}_{x_n}(2,:)\\\mathbf{F}_{u_1}(2,:)\\\mathbf{F}_{u_2}(2,:)\\\vdots\\\mathbf{F}_{u_n}(2,:)\end{array}\right],
% $$
% 
% leading to an $n+m$ by $r$ matrix. Note that due to the norm-1,
% specifying both columns of the monomial factor matrices is redundant.
%
% <html><h3>Phimat - Parameter matrix</h3></html> 
%
% Parameter matrix, which in the context of the MTI, contains the parameters
% for each rank element of the tensor. For example for the state transition tensor
% parameter $\mathrm{F}$ for explicit MTI systems, it would lead to a $n$
% by $r$ factor matrix $\mathbf{F}_{\phi}$.
%
% <html><h2>Output Arguments</h2></html>
%
% <html><h3>CPN1tens - Output tensor</h3></html> 
%
% Output tensor returned as:
%
% * A CP norm-1 (|CPN1|) tensor object.
%
%% Properties
% 
% <html><h3>U - Structure matrix</h3></html> 
% 
% Structure matrix, which for MTI parameter tensors correspond to the
% monomial dimension by the rank $r$.
%
% <html><h3>phi - Parameter matrix</h3></html> 
% 
% Parameter matrix, which for MTI parameter tensors correspond to the
% parameters of each equation by the rank $r$.
%
%% Examples
%
% <html><h3>SISO Multilinear Time-Invariant State-Space Model</h3></html>
%
% Create the state-transition tensor of a SISO MTI state-space model defined by the following system
% equations:
%
% $$
% \dot{x_1} = x_1 x_2 + 0.5 u
% $$
%
% $$
% \dot{x_2} = 2 x_1 u + 7
% $$
%
% The model has the following parameter matrix
%
% $$
% \mathbf{F} = \left(\begin{array}{cccccccc}
% 0 & 0 & 0 & 1 & 0.5 & 0 & 0 & 0 \\ 
% 7 & 0 & 0 & 0 & 0 & 2 & 0 & 0
% \end{array}\right),
% $$
%
% and monomial vector
% 
% $$
% \mathbf{m}\left(x_1,x_2,u\right) = \left(\begin{array}{c}
% 1 \\ x_1 \\ x_2 \\ x_1 x_2 \\ u \\ x_1 u \\ x_2 u \\ x_1 x_2 u
% \end{array}\right).
% $$
%
% The same parameter matrix can be presented in |CPN1| format as:
%
% $$
% \mathbf{F}_\mathrm{U} = \left(\begin{array}{cccc}
% 1 & 0 & 1 & 0 \\ 
% 1 & 0 & 0 & 0 \\
% 0 & 1 & 1 & 0
% \end{array}\right),
% $$
%
% $$
% \mathbf{F}_\mathrm{phi} = \left(\begin{array}{cccc}
% 1 & 0.5 & 0 & 0 \\ 
% 0 & 0 & 2 & 7
% \end{array}\right),
% $$
%
% Initializing the $\mathrm{F}$ tensor in |CPN1| format
F_U = [1  0  1  0;
       1  0  0  0;
       0  1  1  0];
F_phi = [1  0.5 0  0;
         0  0   2  7];
F=CPN1(F_U,F_phi);

%% 
% <html><h3>Normalizing CP tensor to CPN1 tensor object</h3></html>
%
% Create a ktensor for the mti system above with factor matrices
Fx1=[0 1 0 1;1 0 1 0];
Fx2=[0 1 1 1;1 0 0 0];
Fu=[1 0 0 1;0 1 1 0];
FPhi=[1 0.5 0 0;0 0 2 7];
lambda=[1;1;1;1];
CPtens=ktensor(lambda,Fx1,Fx2,Fu,FPhi);
%% 
%
% Create a CPN1 tensor object by normalization
CPN1tens=CPN1(CPtens);
%% References
%
% [1] Lichtenberg, Gerwald; Pangalos, Georg; Cateriano Yáñez, Carlos; Luxa,
% Aline; Jöres, Niklas; Schnelle, Leona; Kaufmann, Christoph (2022): 
% Implicit multilinear modeling. In at - Automatisierungstechnik 70 (1), 
% pp. 13–30. DOI: 10.1515/auto-2021-0133.
%
% <html>
% <style>
% a.button {
%   border-radius: 5px;
%   border: none;
%   color: white;
%   padding: 15px 32px;
%   text-align: center;
%   text-decoration: none;
%   display: inline-block;
%   font-size: 16px;
%   margin: 2px 1px;
%   cursor: pointer;
%   background-color: #0076A8; /* Matlab Blue */
%  }
% a.button:hover {
%  background-color: #3391B9;
% }
% </style>
% <a href="mssDoc.html" class="button">Open mss</a>
% </html>