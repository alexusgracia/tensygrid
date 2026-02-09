function sys = xelim(sys, elim)
%XELIM Eliminates states with the index elim by the matched DC gain method.
%   For continuous-time systems the state derivative are set to zero
%   and the corresponding states are transfered into algebraic variables.
%   For discrete-time systems the state prediction and current state is
%   replaced by the same algebraic variable.

arguments
    sys dmss
    elim (:,1) {mustBeInteger}
end

if sys.ts == 0
    % Continuous-time system
    columnsToEliminate = [];
    for k = 1:length(elim)
        [~,col_t] = find(sys.H.F.stateDerivative.t(elim(k),:));
        [~,col_f] = find(sys.H.F.stateDerivative.f(elim(k),:));
        [~,col_c] = find(sys.H.F.stateDerivative.c(elim(k),:));
    
        columnsToEliminate = [columnsToEliminate; col_t; col_f; col_c];
    end
    
    columnsToEliminate = unique(columnsToEliminate);
    
    sys = eliminateColumns(sys, columnsToEliminate);
    sys = convertStatesToAlgebraics(sys, elim);
else
    % Discrete-time system
    k=1;
    while k <= length(elim)
        [~,col_t] = find(sys.H.F.stateDerivative.t(elim(k),:));
        [~,col_f] = find(sys.H.F.stateDerivative.f(elim(k),:));
        [~,col_c] = find(sys.H.F.stateDerivative.c(elim(k),:));
        colStateDerivative = unique([col_t; col_f; col_c]);
        
        [~,col_t] = find(sys.H.F.state.t(elim(k),:));
        [~,col_f] = find(sys.H.F.state.f(elim(k),:));
        [~,col_c] = find(sys.H.F.state.c(elim(k),:));
        colState = unique([col_t; col_f; col_c]);

        if ismember(colStateDerivative, colState)
            warning(fprintf("State %d cannot be eliminated, since it is multiplied by its predicition (xp%d * x%d) in some term.", elim(k), elim(k), elim(k)))
            elim(k) = [];
        else
            sys = addStateToAlgebraics(sys, elim(k));
            k = k+1;
        end
    end

    sys = eliminateRemainingStates(sys, elim);

end

sys = sys.trivialReduction();


    function sys = eliminateColumns(sys, columnsToEliminate)
        sys.H.F.stateDerivative.t(:, columnsToEliminate) = [];
        sys.H.F.stateDerivative.f(:, columnsToEliminate) = [];
        sys.H.F.stateDerivative.c(:, columnsToEliminate) = [];

        sys.H.F.state.t(:, columnsToEliminate) = [];
        sys.H.F.state.f(:, columnsToEliminate) = [];
        sys.H.F.state.c(:, columnsToEliminate) = [];

        sys.H.F.input.t(:, columnsToEliminate) = [];
        sys.H.F.input.f(:, columnsToEliminate) = [];
        sys.H.F.input.c(:, columnsToEliminate) = [];

        sys.H.F.algebraic.t(:, columnsToEliminate) = [];
        sys.H.F.algebraic.f(:, columnsToEliminate) = [];
        sys.H.F.algebraic.c(:, columnsToEliminate) = [];

        sys.H.F.boolean.t(:, columnsToEliminate) = [];
        sys.H.F.boolean.f(:, columnsToEliminate) = [];
        sys.H.F.boolean.c(:, columnsToEliminate) = [];

        if sys.nEq >0
            sys.H.phi.equality.t(:, columnsToEliminate) = [];
            sys.H.phi.equality.f(:, columnsToEliminate) = [];
            sys.H.phi.equality.c(:, columnsToEliminate) = [];
        end

        if sys.nIneq >0
            sys.H.phi.inequality.t(:, columnsToEliminate) = [];
            sys.H.phi.inequality.f(:, columnsToEliminate) = [];
            sys.H.phi.inequality.c(:, columnsToEliminate) = [];
        end
    end

    function sys = convertStatesToAlgebraics(sys, elim)
        sys.H.F.algebraic.t = [sys.H.F.algebraic.t; sys.H.F.state.t(elim,:)];
        sys.H.F.algebraic.f = [sys.H.F.algebraic.f; sys.H.F.state.f(elim,:)];
        sys.H.F.algebraic.c = [sys.H.F.algebraic.c; sys.H.F.state.c(elim,:)];

        sys.H.F.state.t(elim,:) = [];
        sys.H.F.stateDerivative.t(elim,:) = [];
        sys.H.F.state.f(elim,:) = [];
        sys.H.F.stateDerivative.f(elim,:) = [];
        sys.H.F.state.c(elim,:) = [];
        sys.H.F.stateDerivative.c(elim,:) = [];

        sys.p = sys.p + length(elim);
        sys.n = sys.n - length(elim);

        if ~isempty(sys.stateName)
            sys.algebraicName = [sys.algebraicName, sys.stateName(elim)];
            sys.stateName(elim) = [];
        end
        if ~isempty(sys.stateUnit)
            sys.algebraicUnit = [sys.algebraicUnit, sys.stateUnit(elim)];
            sys.stateUnit(elim) = [];
        end
    end

    function sys = addStateToAlgebraics(sys, ix)
        sys.H.F.algebraic.t = [sys.H.F.algebraic.t; sys.H.F.stateDerivative.t(ix,:) + sys.H.F.state.t(ix,:)];
        sys.H.F.algebraic.f = [sys.H.F.algebraic.f; sys.H.F.stateDerivative.f(ix,:) + sys.H.F.state.f(ix,:)];
        sys.H.F.algebraic.c = [sys.H.F.algebraic.c; sys.H.F.stateDerivative.c(ix,:) + sys.H.F.state.c(ix,:)];

        sys.p = sys.p + 1;
        if ~isempty(sys.stateName)
            sys.algebraicName = [sys.algebraicName, sys.stateName(ix)];
        end
        if ~isempty(sys.stateUnit)
            sys.algebraicUnit = [sys.algebraicUnit, sys.stateUnit(ix)];
        end
    end

    function sys = eliminateRemainingStates(sys, elim)
        sys.H.F.state.t(elim,:) = [];
        sys.H.F.stateDerivative.t(elim,:) = [];
        sys.H.F.state.f(elim,:) = [];
        sys.H.F.stateDerivative.f(elim,:) = [];
        sys.H.F.state.c(elim,:) = [];
        sys.H.F.stateDerivative.c(elim,:) = [];

        sys.n = sys.n - length(elim);
        if ~isempty(sys.stateName)
            sys.stateName(elim) = [];
        end
        if ~isempty(sys.stateUnit)
            sys.stateUnit(elim) = [];
        end
    end

end

