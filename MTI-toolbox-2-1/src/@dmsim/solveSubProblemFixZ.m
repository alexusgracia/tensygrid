function [xp, y, G, F, EqXUZ] = solveSubProblemFixZ(sim, solvingOrderSubProblem, EqXUZ, xpInit, yInit, y0, zInit)
%SOLVESUBPROBLEMFIXZ
arguments
    sim dmsim
    solvingOrderSubProblem
    EqXUZ
    xpInit (1,:)
    yInit (1,:)
    y0 (1,:)
    zInit (1,:)
end
% Torben Warnecke - 11/06/2024

xp = xpInit;
y = yInit;

HphiEq = sim.sys.H.phi.equality.c + sim.sys.H.phi.equality.t - sim.sys.H.phi.equality.f;
HphiIneq = sim.sys.H.phi.inequality.c + sim.sys.H.phi.inequality.t - sim.sys.H.phi.inequality.f;
Hphi = [HphiEq; HphiIneq];

HFt = [sim.sys.H.F.stateDerivative.t; sim.sys.H.F.algebraic.t; sim.sys.H.F.boolean.t];
HFf = [sim.sys.H.F.stateDerivative.f; sim.sys.H.F.algebraic.f; sim.sys.H.F.boolean.f];
HFc = [sim.sys.H.F.stateDerivative.c; sim.sys.H.F.algebraic.c; sim.sys.H.F.boolean.c];

% variable indexes of sub problem
varIdx = solvingOrderSubProblem(:,1);
xpIdx = varIdx(varIdx <= sim.sys.n);
yIdx = varIdx(varIdx>sim.sys.n);
yIdx = yIdx(yIdx <= sim.sys.n+sim.sys.p);
zIdx = varIdx(varIdx>sim.sys.n+sim.sys.p);

eqIdx = solvingOrderSubProblem(:,2);
ineqIdx = eqIdx(eqIdx > sim.sys.nEq);
eqIdx = eqIdx(eqIdx <= sim.sys.nEq);

%zi = zInit(zIdx);
for m = 1:length(zIdx)
    n = zIdx(m);
    v = zInit(n-sim.sys.n-sim.sys.p);
    EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
end

if ~isempty([xpIdx; yIdx])
    if solvingOrderSubProblem(1,3) == 0
        % Sub Problem is nonlinear!
        handle = @(x)optiHandle(HFt, HFf, HFc, HphiEq, HphiIneq, [xpIdx; yIdx], eqIdx, EqXUZ, x);

        lb = -inf * ones([1, sim.sys.n+sim.sys.p]);
        ub = inf * ones([1, sim.sys.n+sim.sys.p]);
        lb = [lb([xpIdx; yIdx])];
        ub = [ub([xpIdx; yIdx])];

        guess = [zeros([1,sim.sys.n]), y0];
        guess = [guess([xpIdx; yIdx])];

        options = optimoptions('lsqnonlin','Display','none',...
            'Algorithm', 'levenberg-marquardt',...
            'FiniteDifferenceType', 'central',...
            'SpecifyObjectiveGradient',true,...
            'OptimalityTolerance', eps^(1/2),...%eps^(2/3),... %eps^2,...
            'FunctionTolerance', eps,...%1/10*sqrt(eps),... %eps^2,...
            'StepTolerance', eps^(1/2),...%eps^(2/3),...
            'MaxIterations', max([length(guess)*5e2, 2e3]),...
            'MaxFunctionEvaluations', max([length(guess)*1e3, 4e3]));%,...
        %'StepTolerance', eps);%eps^(1/2));
        %'CheckGradients',true,...
        %'Algorithm', 'levenberg-marquardt',...
        %'ScaleProblem', 'jacobian',...
        %'InitDamping', 1e2,...

        [solution, ~, ~, ~, ~, ~, ~] = lsqnonlin(handle, guess, lb, ub, options);
    elseif solvingOrderSubProblem(1,3) == -1
        % Sub Problem is linear!
            varIdx = [xpIdx; yIdx];
            A = HphiEq(eqIdx,:) * (EqXUZ .* (HFc(varIdx,:) + HFt(varIdx,:) - HFf(varIdx,:))');
            
            F = 1;
            for k = 1:length(varIdx)
                F = F.* (1 - abs(HFc(varIdx(k),:))) .* (1 - HFt(varIdx(k),:));
            end
            b = HphiEq(eqIdx,:) * (EqXUZ .* F');

            solution = A\-b;
    end
    
    [~,solIdx] = ismember(xpIdx, [xpIdx; yIdx]);
    xp(xpIdx) = solution(solIdx);
    [~,solIdx] = ismember(yIdx, [xpIdx; yIdx]);
    y(yIdx-sim.sys.n) = solution(solIdx);
    
    varIdx = [xpIdx; yIdx];
    for m = 1:length(varIdx)
        n = varIdx(m);
        v = solution(m);
        EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
    end
end

G = Hphi(ineqIdx,:) * EqXUZ;
F = Hphi(eqIdx,:) * EqXUZ;

% if sum(abs(F)>eps^(1/3))>0
%     disp(output.message)
%     error('Solution inaccurate. Simulation stops.')
% end

    function [Eq, Jacobi] = optiHandle(HFt, HFf, HFc, HphiEq, HphiIneq, idx, eqIdx, EqXUZ, x)
        Hphi = [HphiEq; HphiIneq];

        Jacobi = zeros(length(eqIdx), length(idx));

        EqXUZ0 = EqXUZ;
        EqXUZi = EqXUZ;

        k = 1;
        for o = 1:length(idx)

            EqXUZi = EqXUZi .* (1 + HFt(idx(o),:)'.*HFf(idx(o),:)' + HFt(idx(o),:)'*(x(k)-1) - HFf(idx(o),:)'*x(k));
            EqXUZi = EqXUZi .* (1 - abs(HFc(idx(o),:)') + HFc(idx(o),:)'*x(k));

            % calculate Jacobian of the system to compute gradient of the
            % costfunction:
            % 2-FOLD PRODUCT RULE:
            dHdXi_c = EqXUZ0;
            dHdXi_tf = EqXUZ0;
            q = 1;
            for p = 1:length(idx)
                if p == o
                    dHdXi_c = dHdXi_c .* HFc(idx(p),:)' .* (1 + HFt(idx(p),:)'*(x(q)-1) - HFf(idx(p),:)'*x(q) + HFt(idx(p),:)'.*HFf(idx(p),:)');
                    dHdXi_tf = dHdXi_tf .* (HFt(idx(p),:)'-HFf(idx(p),:)') .* (1 - abs(HFc(idx(p),:)') + HFc(idx(p),:)'*x(q));
                    q = q+1;
                else
                    dHdXi_c = dHdXi_c .* (1 - abs(HFc(idx(p),:)') + HFc(idx(p),:)'*x(q));
                    dHdXi_tf = dHdXi_tf .* (1 + HFt(idx(p),:)'.*HFf(idx(p),:)' + HFt(idx(p),:)'*(x(q)-1) - HFf(idx(p),:)'*x(q));
                    q = q+1;
                end
            end
            dHdXi = dHdXi_c + dHdXi_tf;
            Jacobi(:,o) = Hphi(eqIdx,:) * dHdXi;

            k = k+1;
        end
        Eq = Hphi(eqIdx,:) * EqXUZi;
    end
end

