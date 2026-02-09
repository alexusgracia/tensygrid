function [symEq] = symbolicEquations(sys)
%SYMBOLICEQUATIONS Calculate symbolic equations of a dmss model with
% symbolic variables

% Torben Warnecke - 11/06/2024

xp = sym('xp', [sys.n,1]);
x = sym('x', [sys.n,1]);
y = sym('y', [sys.p,1]);
u = sym('u', [sys.m,1]);
z = sym('z', [sys.q,1]);

[Eq, Ineq, ~] = sys.createEquations(xp,x,y,u,z);
symEq = [Eq == 0; 0 >= Ineq];
%symEq = [symEq; Con == 0];
end

