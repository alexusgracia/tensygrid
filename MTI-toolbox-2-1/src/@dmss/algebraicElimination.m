function sys = algebraicElimination(sys, analysisPoints)
%ALGEBRAICREDUCTION Eliminate algrbaic equalities and variables of a dmss model, if its possible
% while maintaining the implicit multilinear structure.
% Indices in analysisPoints won't be reduced, e.g. if they are outputs of
% the model.
%
% For detailed documentation see <a href="matlab:open((which('algebraicEliminationDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024


if nargin < 2
    analysisPoints = [];
end

sys.H = sys.H.convert2doubles;
% sys = sys.normalizeParameters(1e0);

I = sys.incidenceMatrix;
Ialgebraic = I.equality.algebraic;

phi = [sys.H.phi.equality.c + sys.H.phi.equality.t - sys.H.phi.equality.f;...
    sys.H.phi.inequality.c + sys.H.phi.inequality.t - sys.H.phi.inequality.f];

F = [sys.H.F.stateDerivative.c + sys.H.F.stateDerivative.t - 0.5.*sys.H.F.stateDerivative.f;...
    sys.H.F.state.c + sys.H.F.state.t - 0.5.*sys.H.F.state.f;...
    sys.H.F.input.c + sys.H.F.input.t - 0.5.*sys.H.F.input.f;...
    sys.H.F.algebraic.c + sys.H.F.algebraic.t - 0.5.*sys.H.F.algebraic.f;...
    sys.H.F.boolean.c + sys.H.F.boolean.t - 0.5.*sys.H.F.boolean.f];

reducedVar = [];
reducedEq = [];

% go through all algebraic variables (startiung from the last one)
n = 0;
while n < sys.p
    varIdx = sys.p-n;

    % check if the variable is an analysis point
    if ~ismember(varIdx, analysisPoints)
        % if variable is no analysis point, check in which equation it is
        eq = find(Ialgebraic(:,varIdx), sys.nEq);
        eqIdx = eq(end); 
        % Take the last equation, and check if it can be reduced (staying in a multilinear form)

        % check reducability
        isReducable = checkReducability(F, phi, varIdx, eqIdx);

        if isReducable == 1
            % Reduction method (algebraic and trivial Reduction combined)
            [F, phi] = reduction(F, phi, varIdx, eqIdx);
            reducedVar = [reducedVar, varIdx];
            reducedEq = [reducedEq, eqIdx];

            sys.nEq = sys.nEq - 1;
            sys.p = sys.p-1;

            % compute new incidence matrix of the reduced system after
            % each reduction
            Ialgebraic = spones(spones(phi) * spones(F(2*sys.n+sys.m+1:2*sys.n+sys.m+sys.p,:))');
        else
            n = n+1;
        end
    else
        n = n+1;
    end
end

%check if there are all zeros columns in phi
zeroCol = find(all(phi==0));

phi(:,zeroCol) = []; % delete all zeros columns 

sys.H.phi.equality = phi(1:sys.nEq,:);
sys.H.phi.inequality = phi(sys.nEq+1:end,:);

F(:,zeroCol) = []; % delete all zeros columns

sys.H.F.stateDerivative = F(1:sys.n,:);
sys.H.F.state = F(sys.n+1:2*sys.n,:);
sys.H.F.input = F(2*sys.n+1:2*sys.n+sys.m,:);
sys.H.F.algebraic = F(2*sys.n+sys.m+1:2*sys.n+sys.m+sys.p,:);
sys.H.F.boolean = F(2*sys.n+sys.m+sys.p+1:end,:);

if ~isempty(sys.algebraicName)
    sys.algebraicName(reducedVar) =[];
end
if ~isempty(sys.algebraicUnit)
    sys.algebraicUnit(reducedVar) =[];
end

fprintf('System reduced algebraicly by %d algebraic variable(s) and %d equation(s).\n', length(reducedVar), length(reducedEq))

%sys = sys.trivialReduction();

    function isReducable = checkReducability(F, phi, varIdx, eqIdx)
        phiOthers = phi;
        phiOthers(eqIdx,:) = [];
        phiEq = phi(eqIdx, :);

        Fvar = F(2*sys.n+sys.m+varIdx, :);
        FOthers = spones(F);
        FOthers(2*sys.n+sys.m+varIdx,:) = [];
        
        reducable = 0;
        k = 1;
        while (k <= size(phi,2)) && all(reducable==0)
            m = 1;
            while (m <= size(phi,2)) && all(reducable==0)
                reducable = reducable + phiEq(m) * phiOthers(:,k) * ((1-abs(Fvar(k))) .* Fvar(m) - Fvar(k) .* (1-abs(Fvar(m)))) * (FOthers(:,k)' * FOthers(:,m));
%                 if any(reducable ~= 0)
%                     disp("!")
%                 end
                m = m+1;
            end
            k = k+1;
        end

        isReducable = all(reducable==0);
    end

    function [FNew, phiNew] = reduction(F, phi, varIdx, eqIdx)
        phiEq = phi(eqIdx,:);
        %phiEq = find(sign(phiEq),1,'first') * phiEq;

        FVar = F(2*sys.n+sys.m+varIdx,:);

        idxNonzero = find(phiEq);
        idx = find(FVar(idxNonzero), 1);
        phiEq = sign(phiEq(idxNonzero(idx))) * phiEq;

        phiBar = phi;
        phiBar(eqIdx,:) = [];

        FBar = F;
        FBar(2*sys.n+sys.m+varIdx,:) = [];

        FNew = nan(size(FBar,1),1);
        phiNew = [];

        % Algebraic and Trivial Reduction in one!
        q = 1;
        for k = 1:size(phi,2)
            for m = 1:size(phi,2)
                phiTilde = phiBar(:,k) * phiEq(1,m) * ( (1-abs(FVar(k)))*FVar(m) - FVar(k)*(1-abs(FVar(m))) );

                if any(phiTilde~=0)
                    FTilde = FBar(:,k) + FBar(:,m);
                    [Lia, Locb] = ismember(FNew',FTilde','rows'); %check if FTilde already exist in FNew (to not produce duplicates)
                    if ~Lia
                        % if FTilde is new to FNew, the Rank of the
                        % CP-decomposition will be increased
                        phiNew(:,q) = phiTilde;
                        FNew(:,q) = FTilde;
                        q = q+1;
                    else
                        % if FTilde is already in FNew, the according
                        % phiNew column will be extended
                        idx = find(Locb);
                        phiNew(:,idx) = phiNew(:,idx) + phiTilde;
                    end

                end
            end
        end
        phiNew = phiNew ./ max(abs(phiNew),[], 2);
    end
end

