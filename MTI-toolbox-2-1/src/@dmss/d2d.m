function sys = d2d(sys, Ts)
% D2D - Resample discrete-time dynamic iMTI system.
%
%   SYSD2 = d2d(SYSD1,TS) computes a resampled discrete-time model with 
%   using one step forward euler approximation
%
% For detailed documentation see <a href="matlab:open((which('d2dDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024


% Resample discrete-time dmss model
if sys.ts == Ts
    error('dmss model already has a sampling time of %f !', Ts)
end
if isa(sys.ts,'double')
    if sys.ts == 0
        error('dmss model is in continous-time!')
    end
end

%Number of Columns of H.F
N = size(sys.H.F.stateDerivative.c, 2);
H = sys.H;
TsOld = sys.ts;
sys.ts = Ts;

% H.F row of the new auxcilary variable replacing current xp / dot(x)
H.F.algebraic.c = [H.F.algebraic.c, zeros(sys.p, 3*sys.n);...
    H.F.stateDerivative.c, zeros(sys.n, 3*sys.n)];
H.F.algebraic.t = [H.F.algebraic.t, zeros(sys.p, 3*sys.n);...
    H.F.stateDerivative.t, kron(eye(sys.n), [0 0 1])];
H.F.algebraic.f = [H.F.algebraic.f, zeros(sys.p, 3*sys.n);...
    H.F.stateDerivative.f, zeros(sys.n, 3*sys.n)];

% replace H.F row of state Derivative by zeros
H.F.stateDerivative.c = [zeros(sys.n, N + 3*sys.n)];
H.F.stateDerivative.t = [zeros(sys.n, N), kron(eye(sys.n), [1 0 0])];
H.F.stateDerivative.f = [zeros(sys.n, N + 3*sys.n)];

% add additional columns to H.F of the states
H.F.state.c = [H.F.state.c, zeros(sys.n, 3*sys.n)];
H.F.state.t = [H.F.state.t, kron(eye(sys.n), [0 1 0])];
H.F.state.f = [H.F.state.f, zeros(sys.n, 3*sys.n)];

% add addtional columns to H.F of boolean variables
H.F.boolean.c = [H.F.boolean.c, zeros(sys.q, 3*sys.n)];
H.F.boolean.t = [H.F.boolean.t, zeros(sys.q, 3*sys.n)];
H.F.boolean.f = [H.F.boolean.f, zeros(sys.q, 3*sys.n)];

% add addtional columns to H.F of input variables
H.F.input.c = [H.F.input.c, zeros(sys.m, 3*sys.n)];
H.F.input.t = [H.F.input.t, zeros(sys.m, 3*sys.n)];
H.F.input.f = [H.F.input.f, zeros(sys.m, 3*sys.n)];

% add additional equations (row and columns) to H.phi (using forward Euler)
H.phi.equality.c = [H.phi.equality.c, zeros(sys.nEq, 3*sys.n);...
                    zeros(sys.n, N), kron(eye(sys.n), [TsOld Ts-TsOld -Ts])];
H.phi.equality.t = [H.phi.equality.t, zeros(sys.nEq, 3*sys.n);...
                    zeros(sys.n, N), kron(eye(sys.n), [0 0 0])];
H.phi.equality.f = [H.phi.equality.f, zeros(sys.nEq, 3*sys.n);...
                    zeros(sys.n, N), kron(eye(sys.n), [0 0 0])];

H.phi.inequality.c = [H.phi.inequality.c, zeros(sys.nIneq, 3*sys.n)];
H.phi.inequality.t = [H.phi.inequality.t, zeros(sys.nIneq, 3*sys.n)];
H.phi.inequality.f = [H.phi.inequality.f, zeros(sys.nIneq, 3*sys.n)];

sysnew = sys;
sysnew.H = H;
p = sys.p;
sysnew.p = p + sysnew.n;
sysnew.nEq  = sys.nEq + sysnew.n;

if ~isempty(sys.stateName)
    sysnew.stateName = sys.stateName;
end
if ~isempty(sys.algebraicName)
    sysnew.algebraicName = [sys.algebraicName, 'auxilary ' + sys.stateName];
end
if ~isempty(sys.inputName)
    sysnew.inputName = sys.inputName;
end
if ~isempty(sys.booleanName)
    sysnew.booleanName = sys.booleanName;
end

if ~isempty(sys.stateUnit)
    sysnew.stateUnit = sys.stateUnit;
end
if ~isempty(sys.algebraicUnit)
    sysnew.algebraicUnit = [sys.algebraicUnit, sys.stateUnit];
end
if ~isempty(sys.inputUnit)
    sysnew.inputUnit = sys.inputUnit;
end
if ~isempty(sys.booleanUnit)
    sysnew.booleanUnit = sys.booleanUnit;
end

sysnew = sysnew.algebraicElimination(1:p);
%sysnew = sysnew.normalizeParameters(1e0); % is already in algebraicEliminiation

sys = sysnew;
end