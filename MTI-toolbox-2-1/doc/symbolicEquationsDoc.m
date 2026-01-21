%% smybolicEquations
% Outputs the mss-model as symbolic equations
%% Syntax 
%   symbolicEquations(msys)
%   [xp,y]=symbolicEquations(msys)
%% Description
% |symbolicEquations(msys)| Outputs the state equations of the multilinear state-space object
% |msys| as a symbolic equation.
%
% |[xp,y]=symbolicEquations(msys)| Outputs the state equations |xp| and the output equations |y| of the multilinear state-space object
% |msys| as a symbolic equation.
%
% Note that this command requires the Symbolic Math Toolbox from Matlab. 

%% Input Arguments
% |msys| multilinear state-space model of type mss in CPN1 format
%

%% Output Arguments
% |xp| vector or scalar of the state equations as symbolic equations 
%
% |y| vector or scalar of the output equations as symbolic equations 
%% Example 1   
%  Firstly we can create a random continuous-time system with n=2 states, m=1 input, p=2 output, and
%  a rank of r=4:
   msys=rmss(2,1,2,4)
%%   
% Then we can check the equations by:
%
   [xp,y]=symbolicEquations(msys)

%% Example 2
% Assume we have the following second-order explicit MTI model with 2 states and 1 input $$\dot{\mathbf{x}}=\left(\begin{array}{cc} x_2 + x_1x_2\\ 3u_1 + 2x_2 - 3u_1x_1 + 2x_1x_2
% \end{array}\right)$$. We write the explicit MTI model as a |CPN1| object by defining the structural matrix S to form the monomial 
%
    S=[0.5 -0.5; 1 0; 0 1];
%%
% and the parameter matrix with corresponding coefficients of the summands:
%
    phi=[2 0; 4 6];
%%
% Then we can create explicit MTI model as a |CPN1| object and the
% |mss|-object:
%
    obj=CPN1(S,phi);
    msys=mss(obj);
%% 
% To display the equations, simply use:
%
    [xp]=symbolicEquations(msys)
%%
% To get the expanded form, we can use the |expand()| command of the
% Symbolic Math Toolbox.
% 
    expand(xp)

%% See also 
% <mssDoc.html mss>, <rmssDoc.html rmss>