function sysOut = append(sys, varargin)
% Append - Append multiple dmss-objects (diagonally)
%
%   Examples:
%   sys = sys1.append(sys2, sys3, sys4);
%
%   sys = append(sys1, sys2, sys3, sys4);
%
% For detailed documentation see <a href="matlab:open((which('appendDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

R = max([sys.H.sparseIndices.equality.nCols, sys.H.sparseIndices.inequality.nCols]);
nEq = sys.nEq;
nIneq = sys.nIneq;
n = sys.n;
m = sys.m;
p = sys.p;
q = sys.q;

for k = 1:nargin-1
    sys2 = varargin{1,k};
    if class(sys2) ~= "dmss"
        error('the systems attempting to append must be a dmss class object!')
    end

    if class(sys) ~= class(sys2)
        error('The systems attempting to append must be from the same class!')
    end

    % check for similiar steptime!
    if sys.ts ~= sys2.ts
        error('The systems attempting to append must have the same time step size!\n(Reconfigure the systems or try to use c2d/d2c/d2d)')
    end

    R = R + max([sys2.H.sparseIndices.equality.nCols, sys2.H.sparseIndices.inequality.nCols]);
    nEq = nEq + sys2.nEq;
    nIneq = nIneq + sys2.nIneq;
    n = n + sys2.n;
    m = m + sys2.m;
    p = p + sys2.p;
    q = q + sys2.q;
end

sysOut = sys;
sysOut.n = n;
sysOut.m = m;
sysOut.p = p;
sysOut.q = q;
sysOut.nEq = nEq;
sysOut.nIneq = nIneq;

sysOut.stateName = strings(n,1);
sysOut.stateUnit = strings(n,1);
sysOut.inputName = strings(m,1);
sysOut.inputUnit = strings(m,1);
sysOut.algebraicName = strings(p,1);
sysOut.algebraicUnit = strings(p,1);
sysOut.booleanName = strings(q,1);
sysOut.booleanUnit = strings(q,1);

sysOut.manipulatedVariable = [];
sysOut.measuredDisturbance = [];
sysOut.unmeasuredDisturbance = [];

sysOut.H.F.stateDerivative = preallocateTrueFalseMatrix(n, R);
sysOut.H.F.state = preallocateTrueFalseMatrix(n, R);
sysOut.H.F.input = preallocateTrueFalseMatrix(m, R);
sysOut.H.F.algebraic = preallocateTrueFalseMatrix(p, R);
sysOut.H.F.boolean = preallocateTrueFalseMatrix(q, R);

sysOut.H.phi.equality = preallocateTrueFalseMatrix(nEq, R);
sysOut.H.phi.inequality = preallocateTrueFalseMatrix(nIneq, R);

Ri = 0;
nEqi = 0;
nIneqi = 0;
ni = 0;
mi = 0;
pi = 0;
qi = 0;
for k = 0:nargin-1
    if k==0
        sys2 = sys;
    else
        sys2 = varargin{1,k};
    end
    R2 = max([sys2.H.sparseIndices.equality.nCols, sys2.H.sparseIndices.inequality.nCols]);

    sysOut.H.F.stateDerivative = insertDiagonally(sysOut.H.F.stateDerivative, sys2.H.F.stateDerivative, ni+1:ni+sys2.n, Ri+1:Ri+R2);
    sysOut.H.F.state = insertDiagonally(sysOut.H.F.state, sys2.H.F.state, ni+1:ni+sys2.n, Ri+1:Ri+R2);
    sysOut.H.F.input = insertDiagonally(sysOut.H.F.input, sys2.H.F.input, mi+1:mi+sys2.m, Ri+1:Ri+R2);
    sysOut.H.F.algebraic = insertDiagonally(sysOut.H.F.algebraic, sys2.H.F.algebraic, pi+1:pi+sys2.p, Ri+1:Ri+R2);
    sysOut.H.F.boolean = insertDiagonally(sysOut.H.F.boolean, sys2.H.F.boolean, qi+1:qi+sys2.q, Ri+1:Ri+R2);

    sysOut.H.phi.equality = insertDiagonally(sysOut.H.phi.equality, sys2.H.phi.equality, nEqi+1:nEqi+sys2.nEq, Ri+1:Ri+R2);
    sysOut.H.phi.inequality = insertDiagonally(sysOut.H.phi.inequality, sys2.H.phi.inequality, nIneqi+1:nIneqi+sys2.nIneq, Ri+1:Ri+R2);

    if length(sys2.stateName) == sys2.n
        sysOut.stateName(ni+1:ni+sys2.n) = reshape(sys2.stateName, 1, []);
    elseif ~isempty(sys2.stateName)
        error("Missing state names in system %d.", k+1)
    end
    if ~isempty(sys2.stateUnit)
        sysOut.stateUnit(ni+1:ni+sys2.n) = reshape(sys2.stateUnit, 1, []);
    end
    if  length(sys2.inputName) == sys2.m
        sysOut.inputName(mi+1:mi+sys2.m) = reshape(sys2.inputName, 1, []);
    elseif ~isempty(sys2.inputName)
        error("Missing input names in system %d.", k+1)
    end
    if ~isempty(sys2.inputUnit)
        sysOut.inputUnit(mi+1:mi+sys2.m) = reshape(sys2.inputUnit, 1, []);
    end
    if length(sys2.algebraicName) == sys2.p
        sysOut.algebraicName(pi+1:pi+sys2.p) = reshape(sys2.algebraicName, 1, []);
    elseif ~isempty(sys2.algebraicName)
        error("Missing algebraic variable names in system %d.", k+1)
    end
    if ~isempty(sys2.algebraicUnit)
        sysOut.algebraicUnit(pi+1:pi+sys2.p) = reshape(sys2.algebraicUnit, 1, []);
    end
    if length(sys2.booleanName) == sys2.q
        sysOut.booleanName(qi+1:qi+sys2.q) = reshape(sys2.booleanName, 1, []);
    elseif ~isempty(sys2.booleanName)
        error("Missing boolean variable names in system %d.", k+1)
    end
    if ~isempty(sys2.booleanUnit)
        sysOut.booleanUnit(qi+1:qi+sys2.q) = reshape(sys2.booleanUnit, 1, []);
    end

    if ~isempty(sys2.manipulatedVariable)
        sysOut.manipulatedVariable = [sysOut.manipulatedVariable, reshape(sys2.manipulatedVariable+ni, 1, [])];
    end
    if ~isempty(sys2.measuredDisturbance)
        sysOut.measuredDisturbance = [sysOut.measuredDisturbance, reshape(sys2.measuredDisturbance+ni, 1, [])];
    end
    if ~isempty(sys2.unmeasuredDisturbance)
        sysOut.unmeasuredDisturbance = [sysOut.unmeasuredDisturbance, reshape(sys2.unmeasuredDisturbance+ni, 1, [])];
    end

    ni = ni+sys2.n;
    mi = mi+sys2.m;
    pi = pi+sys2.p;
    qi = qi+sys2.q;
    nEqi = nEqi+sys2.nEq;
    nIneqi = nIneqi+sys2.nIneq;
    Ri = Ri+R2;
    
end

    function Mout = preallocateTrueFalseMatrix(nRow, nCol)
        Mout.t = sparse(zeros(nRow, nCol));
        Mout.f = sparse(zeros(nRow, nCol));
        Mout.c = sparse(zeros(nRow, nCol));        
    end
    function Mout = insertDiagonally(Mout, M2, idxRows, idxCols)
        Mout.t(idxRows, idxCols) = M2.t;
        Mout.f(idxRows, idxCols) = M2.f;
        if isempty(M2.c)
            
        elseif isnumeric(M2.c)
            Mout.c(idxRows, idxCols) = M2.c;
        else
            Mout.c=sym(Mout.c);
            Mout.c(idxRows, idxCols) = M2.c;
        end
    end
end

