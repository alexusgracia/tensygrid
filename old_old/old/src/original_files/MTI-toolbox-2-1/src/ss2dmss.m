function dmsys = ss2dmss(lsys)
%SS2DMSS
% Convert linear time invariant state-space model into implicit multilinear
% time invariant state-space model.

% arguments
%  lsys 
% end

if isa(lsys, 'ss')
    E = eye(size(lsys.A,1));
elseif isa(lsys, 'dss')
    E = lsys.E;
else
    error("Only linear state space models stored as ss or dss object are convertible.")
end
A = lsys.A;
B = lsys.B;
C = lsys.C;
D = lsys.D;

n = size(A,1);
m = max(size(B,2), size(D,2));
p = max(size(C,1), size(D,1));

R = 2*n+m+p;

H = hyCPN1();

H.F.stateDerivative = [eye(n), zeros(n,R-n)];
H.F.state = [zeros(n,n), eye(n) ,zeros(n,R-2*n)];
H.F.input = [zeros(m, 2*n), eye(m), zeros(m,R-2*n-m)];
H.F.algebraic = [zeros(p, 2*n+m), eye(p)];
H.F.boolean = double.empty(0, R);

PhiEq = [E, -A, -B, zeros(n,p);...
                zeros(p,n), -C, -D, eye(p)];
H.phi.equality = PhiEq./max(abs(PhiEq),[],2)
H.phi.inequality = double.empty(0, R);

dmsys = dmss(H, lsys.Ts);

end

