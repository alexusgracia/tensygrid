function [Eq, Ineq, Con] = createEquations(sys, xp, x, y, u, z)
% Torben Warnecke - 11/06/2024
    function Eqi = hybridStructureSymbolic(Eqi, M, var)
        if class(var) == "symfun"
            var = var(sym('t'));
        end
        if ~isempty(M.t)
            for k = 1:length(var)
                Eqi = Eqi .* (1 + M.t(k,:)'.*M.f(k,:)' + M.t(k,:)'*(var(k)-1) - M.f(k,:)'*var(k));
                Eqi = Eqi .* (1 - abs(M.c(k,:)') + M.c(k,:)'*var(k));
            end
        end
    end
    function Eqi = hybridStructure(Eqi, M, var, isLogical)
        if ~isempty(M.t)
            tRow = M.t(:,1);
            tCol = M.t(:,2);
            if isLogical
                Eqi = Eqi .* accumarray(tCol,var(tRow)',[M.nCols 1],@all,true);
            else
                Eqi = Eqi .* accumarray(tCol,var(tRow)',[M.nCols 1],@prod,1);
            end
        end
        if ~isempty(M.f)
            fRow = M.f(:,1);
            fCol = M.f(:,2);
            if isLogical
                Eqi = Eqi .* accumarray(tCol,1-var(fRow)',[M.nCols 1],@all,true);
            else
                Eqi = Eqi .* accumarray(fCol, 1-var(fRow)',[M.nCols 1],@prod,1);
            end
        end
        if ~isempty(M.c)
            cRow = M.c(:,1);
            cCol = M.c(:,2);
            cData = M.c(:,3);
            cSigns = M.c(:,4);
            Eqi = Eqi .* accumarray(cCol, 1+cData.*(reshape(var(cRow),[],1)+cSigns), [M.nCols 1],@prod,1);
        end
    end

Eqi = 1;
if isnumeric(xp)
    Eqi = hybridStructure(Eqi, sys.H.sparseIndices.stateDerivative, xp, 0);
else
    Eqi = hybridStructureSymbolic(Eqi, sys.H.F.stateDerivative, xp);
end
if isnumeric(x)
    Eqi = hybridStructure(Eqi, sys.H.sparseIndices.state, x, 0);
else
    Eqi = hybridStructureSymbolic(Eqi, sys.H.F.state, x);
end
if isnumeric(y)
    Eqi = hybridStructure(Eqi, sys.H.sparseIndices.algebraic, y, 0);
else
    Eqi = hybridStructureSymbolic(Eqi, sys.H.F.algebraic, y);
end
if isnumeric(u)
    Eqi = hybridStructure(Eqi, sys.H.sparseIndices.input, u, 0);
else
    Eqi = hybridStructureSymbolic(Eqi, sys.H.F.input, u);
end
if isnumeric(z)
    Eqi = hybridStructure(Eqi, sys.H.sparseIndices.boolean, logical(z), 1);
else
    Eqi = hybridStructureSymbolic(Eqi, sys.H.F.boolean, z);
end

    function Eqn = hybridParameter(Eqi, M)
        if isempty(M.t)
            Eqn = [];
        else
            Eqn = (M.t + M.c - M.f) * Eqi;
        end
    end

Eq = hybridParameter(Eqi, sys.H.phi.equality);
Ineq = hybridParameter(Eqi, sys.H.phi.inequality);
Con = z .* (1-z);
end