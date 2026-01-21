function [xp, y, z, F, G] = solveStepBruteForce(sim, x0, u, solvingOrder, yguess, xpguess, zguess)
arguments
    sim dmsim
    x0 (1,:)
    u (1,:)
    solvingOrder
    yguess (1,:)
    xpguess (1,:)
    zguess  (1,:)
end
% Torben Warnecke - 11/06/2024

xp = NaN(1, sim.sys.n);
y = NaN(1, sim.sys.p);
%z = NaN(1, sim.sys.q);
z = zguess;

EqXUZ = 1;
if sim.sys.n ~=0
    EqXUZ = hybridStructure(EqXUZ, sim.sys.H.sparseIndices.state, x0);
end
if sim.sys.m ~=0
    EqXUZ = hybridStructure(EqXUZ, sim.sys.H.sparseIndices.input, u);
end

HphiEq = sim.sys.H.phi.equality.c + sim.sys.H.phi.equality.t - sim.sys.H.phi.equality.f;
HphiIneq = sim.sys.H.phi.inequality.c + sim.sys.H.phi.inequality.t - sim.sys.H.phi.inequality.f;
Hphi = [HphiEq; HphiIneq];

HFt = [sim.sys.H.F.stateDerivative.t; sim.sys.H.F.algebraic.t; sim.sys.H.F.boolean.t];
HFf = [sim.sys.H.F.stateDerivative.f; sim.sys.H.F.algebraic.f; sim.sys.H.F.boolean.f];
HFc = [sim.sys.H.F.stateDerivative.c; sim.sys.H.F.algebraic.c; sim.sys.H.F.boolean.c];

l = 1;
while (l <= size(solvingOrder,1))
    if solvingOrder(l, 3) == 1
        if solvingOrder(l,1) <= sim.sys.n+sim.sys.p
            n = solvingOrder(l,1);
            e = solvingOrder(l,2);

            v = explicitSolution(Hphi(e,:), EqXUZ, HFc(n,:), HFt(n,:), HFf(n,:));

            EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
            l = l+1;

            if n<=sim.sys.n
                xp(1,n) = v;
            else
                y(1,n - sim.sys.n) = v;
            end
        else
            n = solvingOrder(l,1);
            e = solvingOrder(l,2);

            HphiXU = Hphi(e,:);
            z0 = zguess(n-sim.sys.n-sim.sys.p);
            % z0 = 0;
            EqXUZ0 = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(z0-1) - HFf(n,:)'*z0) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*z0);
            F0 = HphiXU * EqXUZ0;

            if sign(F0) <=0
                v = z0;
            else
                v = ~z0;
            end
            EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
            z(1,n - sim.sys.n - sim.sys.p) = v;
            l = l+1;
        end
    else
        solvingOrderSubProblem_i = solvingOrder(solvingOrder(:,5)==solvingOrder(l,5), :);
        [rowSubProblem, ~] = find(solvingOrder(l:end,5)>solvingOrderSubProblem_i(1,5),1,'first');
        if isempty(rowSubProblem)
            l = size(solvingOrder,1)+1;
        else
            l = l+rowSubProblem-1;
        end

        [xp, y, z, ~, ~, EqXUZ] = sim.solveSubProblemVariableZ(solvingOrderSubProblem_i, EqXUZ, xpguess, y, yguess, z);
    
    end
end

if ~isempty(HphiEq)
    F = HphiEq * EqXUZ;
else
    F = [];
end
if ~isempty(HphiIneq)
    G = HphiIneq * EqXUZ;
else
    G = [];
end

% if any(isnan(F))|| any(isnan(G))
%     %warning('solveStepBruteForce:NoSolution',"No solution could be found.")
% end
% if any(abs(F)>eps^(1/4))
%     warning('solveStepBruteForce:Inaccurate',"Simulation might be inaccurate. Absolute function value is bigger then eps^(1/4)." + sprintf(' Maximum absolute function value is %f\n',max(abs(F))))
% end

%% Functions
    function Eqi = hybridStructure(Eqi, M, var)
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

    function v = explicitSolution(HphiXU, EqXUZ, HFcn, HFtn, HFfn)
        B = HphiXU * ((1-abs(HFcn + HFtn.*HFfn - HFtn))' .* EqXUZ);
        A = HphiXU * ((HFtn - HFfn + HFcn)' .* EqXUZ);
        v = -B/A;
    end
end

