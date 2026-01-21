function [H, ts] = mss2dmss(sys)
%MSS2DMSS Convert mss-object into dmss-object
%   converts an explicit multilinear timeinvariant model into an implicit
%   multilinear timeinvariant model

% Torben Warnecke - 11/06/2024

arguments
    sys mss 
end

H = hyCPN1();
N = size(sys.F.U,2);
M = size(sys.G.U,2);
n = sys.n;
p = sys.p;
m = sys.m;

H.F.stateDerivative = [diag(ones([n,1])), zeros([n, p+N+M])];
H.F.algebraic = [zeros([p,n]), diag(ones([p,1])), zeros([p, N+M])];
H.F.state = [zeros([n,n+p]), sys.F.U(1:n,:), sys.G.U(1:n,:)];
H.F.input = [zeros([m,n+p]), sys.F.U(n+1:end,:), sys.G.U(n+1:end,:)];

H.phi.equality = [diag(-ones([n+p,1])), [sys.F.phi, zeros([n,M]); zeros([p,N]), sys.G.phi]];

ts = sys.ts;

end

