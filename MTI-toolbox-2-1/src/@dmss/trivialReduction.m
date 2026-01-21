function sys = trivialReduction(sys)
% REDUCE FACTOR MATRIX
% combine all factor matrices to one big matrices (stacked on top of each
% other (stacked over rows)

% Torben Warnecke - 11/06/2024

if isfield(sys.H.F.state, 't') == 1
    F.c = [sys.H.F.stateDerivative.c; sys.H.F.state.c; sys.H.F.algebraic.c; sys.H.F.input.c; sys.H.F.boolean.c];
    F.t = [sys.H.F.stateDerivative.t; sys.H.F.state.t; sys.H.F.algebraic.t; sys.H.F.input.t; sys.H.F.boolean.t];
    F.f = [sys.H.F.stateDerivative.f; sys.H.F.state.f; sys.H.F.algebraic.f; sys.H.F.input.f; sys.H.F.boolean.f];
    N = size(F.c,1);
    HF = [F.c; F.t; F.f];
else
    HF = [sys.H.F.stateDerivative; sys.H.F.state; sys.H.F.algebraic; sys.H.F.input; sys.H.F.boolean];
end

% find unique columns of all Factor matrices
[Hfactor_transposed, k, l] = unique(HF', 'rows', 'stable');

% rebuild the single factor matrices
F.c = [];
F.c = Hfactor_transposed(:,1:N)';
F.t = [];
F.t = double(Hfactor_transposed(:,N+1:2*N)') == 1;
F.f = [];
F.f = double(Hfactor_transposed(:,2*N+1:3*N)') == 1;

sys.H.F.stateDerivative.c = F.c(1:sys.n,:);
sys.H.F.state.c = F.c(sys.n+1:2*sys.n,:);
sys.H.F.algebraic.c = F.c(2*sys.n+1:2*sys.n+sys.p,:);
sys.H.F.input.c = F.c(2*sys.n+sys.p+1:2*sys.n+sys.p+sys.m,:);
sys.H.F.boolean.c = F.c(2*sys.n+sys.p+sys.m+1:end,:);

sys.H.F.stateDerivative.t = F.t(1:sys.n,:);
sys.H.F.state.t = F.t(sys.n+1:2*sys.n,:);
sys.H.F.algebraic.t = F.t(2*sys.n+1:2*sys.n+sys.p,:);
sys.H.F.input.t = F.t(2*sys.n+sys.p+1:2*sys.n+sys.p+sys.m,:);
sys.H.F.boolean.t = F.t(2*sys.n+sys.p+sys.m+1:end,:);

sys.H.F.stateDerivative.f = F.f(1:sys.n,:);
sys.H.F.state.f = F.f(sys.n+1:2*(sys.n),:);
sys.H.F.algebraic.f = F.f(2*sys.n+1:2*sys.n+sys.p,:);
sys.H.F.input.f = F.f(2*sys.n+sys.p+1:2*sys.n+sys.p+sys.m,:);
sys.H.F.boolean.f = F.f(2*sys.n+sys.p+sys.m+1:end,:);

% REDUCE PARAMETER MATRIX
% combine all parameter matrices
Hphi.equality = sys.H.phi.equality.c + sys.H.phi.equality.t - sys.H.phi.equality.f;
Hphi.inequality = sys.H.phi.inequality.c + sys.H.phi.inequality.t - sys.H.phi.inequality.f;

% combine columns, from nonunique colums of factor matrices
HphiNew.equality = [];
HphiNew.inequality = [];
for n = 1:length(k)
    if ~isempty(Hphi.equality)
        HphiNew.equality = [HphiNew.equality, sum(Hphi.equality(:,l == n),2)];
    end
    if ~isempty(Hphi.inequality)
        HphiNew.inequality = [HphiNew.inequality, sum(Hphi.inequality(:,l == n),2)];
    end
end

% check for duplicats of equations and inequations
[HphiNew.equality, i11, i12] = unique(HphiNew.equality,'rows', 'stable');

[HphiNew.inequality, i21, i22] = unique(HphiNew.inequality,'rows', 'stable');

%fprintf('Reduced by %d column(s) and %d equation(s) with trivial Reduction due to duplications.\n', length(j)-length(i), (length(i12)-length(i11))+(length(i22)-length(i21)))

[Hphi, i31, i32] = unique([HphiNew.equality; HphiNew.inequality],'rows', 'stable');
if length(i31) ~= length(i32)
    disp('At least one inequality is equal to an equation! The inequality has been reduced.')
    N = size(HphiNew.equality,1);
    sys.H.phi.equality = Hphi(1:N,:);
    sys.H.phi.inequality = Hphi(N+1:end,:);
else
    sys.H.phi.equality = HphiNew.equality;
    sys.H.phi.inequality = HphiNew.inequality;
end
if ~isnumeric(sys.H.phi.equality.c)
    try
        sys.H.phi.equality.c = double(sys.H.phi.equality.c);
    catch
    end
end
if ~isnumeric(sys.H.phi.inequality.c)
    try
        sys.H.phi.inequality.c = double(sys.H.phi.inequality.c);
    catch
    end
end

if length(l)-length(k)+length(i12)-length(i11)+length(i22)-length(i21)+length(i32)-length(i31) > 0
    fprintf('Reduced by %d column(s) and %d equation(s) with trivial Reduction due to duplications.\n', length(l)-length(k), length(i12)-length(i11)+length(i22)-length(i21)+length(i32)-length(i31))
end

sys.nIneq = size(sys.H.phi.inequality.c,1);
sys.nEq = size(sys.H.phi.equality.c,1);
end

