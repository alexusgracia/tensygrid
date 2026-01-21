%% MEXaccelerate
% enable MEX acceleration (experimental)
%% Syntax
% MEXaccelerate(mss,enable)
%% Description
% Enable MEX acceleration for simulations. 
%
% Note: Current implementation changes the system in the input. A copy of
% this changed input is returned as output for flexible usage.
%
% For example usage of this feature, see /tests/mexEval.m.
%% Input Arguments 
% |mss|: mss object 
%
% |enable|  : enable flag, logical
%% Example 
A = [-1.5 -3; 3 -1];
B = [1.3; 0];
C = [1.15 2.3];
D = [0];

F = CPN1(diag(ones([3,1])), [A, B]);
G = CPN1(diag(ones([3,1])), [C, D]);
msys = mss(F,G);

x0 = [-0.2 0.3];
t = 0:0.05:8;
u = zeros(length(1),1)';
u(t>=2) = 1;

msys.MEXaccelerate(true);
% msys = MEXaccelerate(msys,true) has the same result as the above syntax


[y, tOut, x] = msim(msys, u', t, x0);

