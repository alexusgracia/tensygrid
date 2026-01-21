function J = jacobian(sys, xp, x, u, y, z)
%JACOBIAN Calculate jacobian matrix of the set of equations of a
%   dmss model.
%
%   The jacobian matrix will be stored as a struct with entries for each
%   variable type. The rows of the sub matrices correspond to the different
%   equations.
%
%   J = jacobian(sys): Calculate jacobian matrix with default symbolic variables [xp, x, u, y].
%
%   J = jacobian(sys, xp, x, u, y): Calculate jacobian matrix for given
%   variable values (numeric or symbolic).
%
%   J = jacobian(sys, xp, x, u, y, z): Calculate jacobian matrix for given
%   variable values (numeric or symbolic) for hybrid dmss model (Experimental).
%
% For detailed documentation see <a href="matlab:open((which('jacobianDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

arguments
    sys dmss
    xp (:, 1) = sym('xp',[sys.n,1])
    x (:, 1) = sym('x',[sys.n,1])
    u (:, 1) = sym('u',[sys.m,1])
    y (:, 1) = sym('y',[sys.p,1])
    z (:, 1) = sym('z',[sys.q,1])
end

if length(xp) ~= sys.n
    if isa(xp, 'dmss')
        error(['When calling a class function, the call can be either "sys.jacobian(xp, x, u, y)" or the object "sys" can be an input argument as in "jacobian(sys, xp, x, u, y)". Both is not possible.' ...
            'Try to delete sys as an input argument: "sys.jacobian(xp, x, u, y)".'])
    else
        error('Number of xp-variables must equal number of state variables (n)!')
    end
elseif length(x) ~= sys.n
    error('Number of x-variables must equal number of state variables (n)!')
elseif length(y) ~= sys.p
    error('Number of y-variables must equal number of algebraic variables (p)!')
elseif length(u) ~= sys.m
    error('Number of u-variables must equal number of input variables (m)!')
elseif length(z) ~= sys.q
    error('Number of z-variables must equal number of boolean variables (q)!')
end

    function Ji = hybridParameter(dHdXi, M)
        if isempty(M.t)
            Ji = [];
        else
            Ji = (M.t + M.c - M.f) * dHdXi;
        end
    end

    function Eqi = hybridStructure(Eqi, M, var, k)
        n = length(var);
        if k>0 && k<=n
            if ~isempty(M.t)
                tRow = M.t(:,1);
                tCol = M.t(:,2);
                Eqi_t = accumarray(tCol, (tRow==k), [M.nCols 1],@sum,0) .* accumarray(tCol, 1+((tRow~=k).*(var(tRow)-1)), [M.nCols 1],@prod,1);
            else
                Eqi_t = 0;
            end
            if ~isempty(M.f)
                fRow = M.f(:,1);
                fCol = M.f(:,2);
                Eqi_f = accumarray(fCol, -(fRow==k), [M.nCols 1],@sum,0) .* accumarray(fCol, (1-(fRow~=k).*var(fRow)), [M.nCols 1],@prod,1);
            else
                Eqi_f = 0;
            end
            if ~isempty(M.c)
                cRow = M.c(:,1);
                cCol = M.c(:,2);
                cData = M.c(:,3);
                cSigns = M.c(:,4);
                Eqi_c = accumarray(cCol, cData.*(cRow==k), [M.nCols 1],@sum,0) .* accumarray(cCol, 1+cData.*(cRow~=k).*(var(cRow)+cSigns), [M.nCols 1],@prod,1);
            else
                Eqi_c = 0;
            end
            Eqi = Eqi .* (Eqi_c+Eqi_t+Eqi_f);
        else
            if ~isempty(M.t)
                tRow = M.t(:,1);
                tCol = M.t(:,2);
                Eqi = Eqi .* accumarray(tCol,var(tRow)',[M.nCols 1],@prod,1);
            end
            if ~isempty(M.f)
                fRow = M.f(:,1);
                fCol = M.f(:,2);
                Eqi = Eqi .* accumarray(fCol, 1-var(fRow)',[M.nCols 1],@prod,1);
            end
            if ~isempty(M.c)
                cRow = M.c(:,1);
                cCol = M.c(:,2);
                cData = M.c(:,3);
                cSigns = M.c(:,4);
                Eqi = Eqi .* accumarray(cCol, 1+cData.*(reshape(var(cRow),[],1)+cSigns), [M.nCols 1],@prod,1);
            end
        end
    end
var = [];
if ~isempty(xp)
    var = [var; xp];
end
if ~isempty(x)
    var = [var; x];
end
if ~isempty(u)
    var = [var; u];
end
if ~isempty(y)
    var = [var; y];
end
% if ~isempty(u)
%     var = [var; u];
% end
if ~isempty(z)
    var = [var; z];
end
%var = [xp; x; y; u; z];

JEqi = [];
JIneqi = [];

for k = 1:length(var)
    if isfield(sys.H.F.state, 't')==1 && isnumeric(var)
        dHdXi = 1;

        dHdXi = hybridStructure(dHdXi, sys.H.sparseIndices.stateDerivative, xp, k);
        dHdXi = hybridStructure(dHdXi, sys.H.sparseIndices.state, x, k-sys.n);
        dHdXi = hybridStructure(dHdXi, sys.H.sparseIndices.input, u, k-2*sys.n);
        dHdXi = hybridStructure(dHdXi, sys.H.sparseIndices.algebraic, y, k-2*sys.n-sys.m);
        dHdXi = hybridStructure(dHdXi, sys.H.sparseIndices.boolean, z, k-2*sys.n-sys.m-sys.p);

    elseif isfield(sys.H.F.state, 't')==1 && ~isnumeric(var)
        F.c = [sys.H.F.stateDerivative.c; sys.H.F.state.c; sys.H.F.input.c; sys.H.F.algebraic.c; sys.H.F.boolean.c];
        F.t = [sys.H.F.stateDerivative.t; sys.H.F.state.t; sys.H.F.input.t; sys.H.F.algebraic.t; sys.H.F.boolean.t];
        F.f = [sys.H.F.stateDerivative.f; sys.H.F.state.f; sys.H.F.input.f; sys.H.F.algebraic.f; sys.H.F.boolean.f];

        % 2-FOLD PRODUCT RULE:
        dHdXi_c = 1;
        dHdXi_tf = 1;
        for l = 1:length(var)
            if l == k
                dHdXi_c = dHdXi_c .* F.c(l,:)' .* (1 + F.t(l,:)'*(var(l)-1) - F.f(l,:)'*var(l) + F.t(l,:)'.*F.f(l,:)');
                dHdXi_tf = dHdXi_tf .* (F.t(l,:)'-F.f(l,:)') .* (1 - abs(F.c(l,:)') + F.c(l,:)'*var(l));
            else
                dHdXi_c = dHdXi_c .* (1 - abs(F.c(l,:)') + F.c(l,:)'*var(l)) .* (1 + F.t(l,:)'*(var(l)-1) - F.f(l,:)'*var(l) + F.t(l,:)'.*F.f(l,:)');
                dHdXi_tf = dHdXi_tf .* (1 - abs(F.c(l,:)') + F.c(l,:)'*var(l)) .* (1 + F.t(l,:)'*(var(l)-1) - F.f(l,:)'*var(l) + F.t(l,:)'.*F.f(l,:)');
            end
        end
        dHdXi = dHdXi_c + dHdXi_tf;

%     else
%         if isnumeric(Data) && isnumeric(var)
%             dHdXi = (accumarray(Col, Data.*(Row==k), [nCols 1],@sum,0) .* accumarray(Col, 1+Data.*(Row~=k).*(var(Row)+Signs), [nCols 1],@prod,1));
%         else
%             dHdXi = 1;
%             for l = 1:length(var)
%                 if l == k
%                     dHdXi = dHdXi .* F(l,:)';
%                 else
%                     dHdXi = dHdXi .* (1 - abs(F(l,:)') + F(l,:)'*var(l));
%                 end
%             end
%         end
    end
    
    if isfield(sys.H.phi.equality, 't')==1
        JEqi = [JEqi, hybridParameter(dHdXi, sys.H.phi.equality)];
        JIneqi = [JIneqi, hybridParameter(dHdXi, sys.H.phi.inequality)];
    else
        JEqi = [JEqi, sys.H.phi.equality * dHdXi];
        JIneqi = [JIneqi, sys.H.phi.inequality * dHdXi];
    end
end

if isempty(JEqi)
    J.equality.stateDerivative = [];
    J.equality.state = [];
    J.equality.algebraic = [];
    J.equality.input = [];
    J.equality.boolean = [];
else
    J.equality.stateDerivative = JEqi(:, 1:sys.n);
    J.equality.state = JEqi(:, sys.n+1:2*sys.n);
    J.equality.input = JEqi(:, 2*sys.n+1:2*sys.n+sys.m);
    J.equality.algebraic = JEqi(:, 2*sys.n+sys.m+1:2*sys.n+sys.m+sys.p);
    J.equality.boolean = JEqi(:, 2*sys.n+sys.p+sys.m+1:2*sys.n+sys.p+sys.m+sys.q);
end

if isempty(JIneqi)
    J.inequality.stateDerivative = [];
    J.inequality.state = [];
    J.inequality.algebraic = [];
    J.inequality.input = [];
    J.inequality.boolean = [];
else
    J.inequality.stateDerivative = JIneqi(:, 1:sys.n);
    J.inequality.state = JIneqi(:, sys.n+1:2*sys.n);
    J.inequality.input = JIneqi(:, 2*sys.n+1:2*sys.n+sys.m);
    J.inequality.algebraic = JIneqi(:, 2*sys.n+sys.m+1:2*sys.n+sys.m+sys.p);
    J.inequality.boolean = JIneqi(:, 2*sys.n+sys.p+sys.m+1:2*sys.n+sys.p+sys.m+sys.q);
end

end

