function [solvingOrder, I, Isorted] = determineParallizedSolvingOrder(sim, withoutZ)
arguments
    sim dmsim
    withoutZ logical
end
%DETERMINESOLVINGORDER

% Torben Warnecke - 11/06/2024

I = sim.sys.incidenceMatrix;
I = [I.equality.stateDerivative, I.equality.algebraic, I.equality.boolean;...
    I.inequality.stateDerivative, I.inequality.algebraic, I.inequality.boolean];
%I = I(:, [1:sim.sys.n, 2*sim.sys.n+sim.sys.m+1:end]);

if withoutZ == 1
    I = I(1:end-sim.sys.nIneq,1:sim.sys.n+sim.sys.p);
end

[m,n] = size(I);
if m~=n 
    error('The system is not a square system! The number of unknown variables do not equal the number of equations and constraints.')
end

[Isorted, sortedRowIdx, sortedColIdx] = sim.sys.orderIncidenceMatrix(I);

n = size(Isorted,1);

%% Identify block-structures (parallizable sub-Sets)
subSet = ones([n,1]);
subProblem = ones([n,1]);
for l = 2:n
    if ~any(Isorted(l:end, 1:l-1), 'all')
        if ~any(Isorted(1:l-1, l:end), 'all')
            subSet(l:end) = subSet(l)+1;
        end
    end
    if ~any(Isorted(1:l-1, l:end), 'all')
        subProblem(l:end) = subProblem(l)+1;
    end
end

% Identify which equations (rows of Jacobian matrix) can be explisitly
% solved for one unknown variable
explicitSolvable = zeros([n,1]);
k = 0;
for l = 1:n
    [~, idx] = find(Isorted(l, l+1:end));
    if (max(idx)+l)>k
        k = max(idx)+l;
    end
    if (l>k)&&(isempty(idx))
        explicitSolvable(l) = 1;
    end
end




solvingOrder = [sortedColIdx, sortedRowIdx, explicitSolvable, subSet, subProblem];

idx = zeros([n,1]);
l=1;
% for k = 2:n
%     if solvingOrder(k,3) ~=1
%         idx(k) = l;
%     elseif solvingOrder(k-1,3)~=1
%         l = l+1;
%     end
% end
implicitSubProblems = unique((solvingOrder(:,3) == 0) .* solvingOrder(:,5));
implicitSubProblems= implicitSubProblems(implicitSubProblems>0);

binaryPhi = [(sim.sys.H.phi.equality.c~=0 + sim.sys.H.phi.equality.t~=0 + sim.sys.H.phi.equality.f~=0)];%;...
            %(sim.sys.H.phi.inequality.c~=0 + sim.sys.H.phi.inequality.t~=0 + sim.sys.H.phi.inequality.f~=0)];
binaryF = [(sim.sys.H.F.stateDerivative.c~=0 + sim.sys.H.F.stateDerivative.t~=0 + sim.sys.H.F.stateDerivative.f~=0);...
            (sim.sys.H.F.algebraic.c~=0 + sim.sys.H.F.algebraic.t~=0 + sim.sys.H.F.algebraic.f~=0)];%;...
            %(sim.sys.H.F.boolean.c~=0 + sim.sys.H.F.boolean.t~=0 + sim.sys.H.F.boolean.f~=0)];

for k = 1:length(implicitSubProblems)
    implicitSetIdx = solvingOrder(:,5)==implicitSubProblems(k);
    
    eqIdx = solvingOrder(implicitSetIdx, 2);
    eqIdx = eqIdx(eqIdx<=sim.sys.nEq);

    varIdx = solvingOrder(implicitSetIdx, 1);
    varIdx = varIdx(varIdx<=sim.sys.n+sim.sys.p);

    columns = sparse(zeros(1, sim.sys.H.R));
    for l =1:length(eqIdx)
        for m = 1:length(varIdx)
            try
                columns = any([columns; binaryPhi(eqIdx(l),:) .* binaryF(varIdx(m),:)]);
            catch
                disp('!')
            end
        end
    end

    numberOfNonzeroEntriesInColumnsOfF = sum(binaryF(varIdx,:),1);
    if all(numberOfNonzeroEntriesInColumnsOfF<=1)
        if all(varIdx <=sim.sys.n+sim.sys.p)
            solvingOrder(implicitSetIdx,3) = -1;
        %else
            %solvingOrder(implicitSetIdx,3) = -2;
        end
    end
end

% SolvingOrder: [UnknownIndex, EquationIndex, explicitMarker, subSet, subProblems]
%   -> explicitMarker: 
%   1 if equation is explicit solvable for unknown, 
%   0 if equations and following equations with 0 need to be solved implicitly
%   and are nonlinear, 
%   -1 if equations and following equations with -1 are
%   implicit but linear equations
end

