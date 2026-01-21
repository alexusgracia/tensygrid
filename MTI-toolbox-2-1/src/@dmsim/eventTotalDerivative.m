function [xp1, y1, z1, dy1, F, G] = eventTotalDerivative(sim, ie, t0, x0, y0, z0, u, ut, solvingOrder, opt) %inequalityTolerance)
% Torben Warnecke - 11/06/2024

% calculate inputs u at the event time step, as well as the derivative of u
% (1-Step Euler)
[~,idx_t] = min(abs(ut - t0));
u0 = zeros(1,sim.sys.m);
du = zeros(1,sim.sys.m);
for l = 1:size(u,2)
    %u0(1,j) = interp1(ut, u(:,j), t0, 'linear');
    if ut(idx_t) ~= t0
        u0(1,l) = interp1q(ut, u(:,l), t0);
        du(1,l) = (u0(1,l) - u(idx_t, l))/(t0 - ut(idx_t));
    elseif idx_t == 1
        u0(1,l) = u(1, l);
        du(1,l) = (u0(1,l) - u(idx_t+1, l))/(t0 - ut(idx_t+1));
    elseif ut(idx_t) == t0
        u0(1,l) = u(idx_t, l);
        du(1,l) = (u0(1,l) - u(idx_t-1, l))/(t0 - ut(idx_t-1));
    end
end


%% Identify the sub-set and its equations and variables
if ~isempty(ie)
% index of the inequality in the solvingOrder-matrix
IneqIdx = ie + sim.sys.nEq;
for k = 1:length(IneqIdx);
    idx(k) = find((solvingOrder(:,2)-IneqIdx(k))==0);
end

% determine the sub-set and sub-problem of concerned equations to reduce complexity
% subSet = min(solvingOrder(idx,4));
% subProblem = min(solvingOrder(idx,5));
subSet = sort(solvingOrder(idx,4));
subProblem = sort(solvingOrder(idx,5));

%subProblemRowIdx = solvingOrder(:,5) == subProblem;
%solvingOrderSubProblem = solvingOrder(subProblemRowIdx,:);
for k = 1:length(subProblem)
    subProblemRowIdx(:,k) = solvingOrder(:,5) == subProblem(k);
    solvingOrderSubProblem{k} = solvingOrder(subProblemRowIdx(:,k),:);
end

%previousSubProblemsRowIdx = 1 : find(solvingOrder(:,5)==subProblem, 1, 'first')-1;
previousSubProblemsRowIdx = 1 : find(solvingOrder(:,5)==min(subProblem), 1, 'first')-1;
solvingOrderPreviousSubProblems = solvingOrder(previousSubProblemsRowIdx,:);

followingSubProblemsRowIdx = find(solvingOrder(:,5)==min(subProblem)+1, 1, 'first') : length(solvingOrder(:,1));
%followingSubProblemsRowIdx = find(solvingOrder(:,5)==subProblem+1, 1, 'first') : length(solvingOrder(:,1));
solvingOrderFollowingSubProblems = solvingOrder(followingSubProblemsRowIdx,:);

else
    solvingOrderPreviousSubProblems = solvingOrder;
    solvingOrderSubProblem = [];
    solvingOrderFollowingSubProblems =[];
end

%% Solve previuous sub problems for fixed (previuous) z values
zInit = z0;

xpInit = nan(1,sim.sys.n);
yInit = nan(1,sim.sys.p);

% xpInit = NaN; %<-
% yInit = NaN; %<-

% xpInit = [];
% yInit = [];

% xpInit = zeros(1,sim.sys.n);
% yInit = zeros(1,sim.sys.p);

EqXUZ = 1;
if sim.sys.n ~=0
    EqXUZ = hybridStructure(EqXUZ, sim.sys.H.sparseIndices.state, x0);
end
if sim.sys.m ~=0
    EqXUZ = hybridStructure(EqXUZ, sim.sys.H.sparseIndices.input, u0);
end

HphiEq = sim.sys.H.phi.equality.c + sim.sys.H.phi.equality.t - sim.sys.H.phi.equality.f;
HphiIneq = sim.sys.H.phi.inequality.c + sim.sys.H.phi.inequality.t - sim.sys.H.phi.inequality.f;
Hphi = [HphiEq; HphiIneq];

HFt = [sim.sys.H.F.stateDerivative.t; sim.sys.H.F.algebraic.t; sim.sys.H.F.boolean.t];
HFf = [sim.sys.H.F.stateDerivative.f; sim.sys.H.F.algebraic.f; sim.sys.H.F.boolean.f];
HFc = [sim.sys.H.F.stateDerivative.c; sim.sys.H.F.algebraic.c; sim.sys.H.F.boolean.c];

l = 1;
while l <= size(solvingOrderPreviousSubProblems,1)
    if solvingOrderPreviousSubProblems(l, 3) == 1
        if solvingOrderPreviousSubProblems(l,1) <= sim.sys.n+sim.sys.p
            n = solvingOrderPreviousSubProblems(l,1);
            e = solvingOrderPreviousSubProblems(l,2);

            v = explicitSolution(Hphi(e,:), EqXUZ, HFc(n,:), HFt(n,:), HFf(n,:));

            EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
            l = l+1;

            if n<=sim.sys.n
                xpInit(1,n) = v;
            else
                yInit(1,n - sim.sys.n) = v;
            end
        else
            n = solvingOrderPreviousSubProblems(l,1);

            v = z0(1, n - sim.sys.n - sim.sys.p);

            EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
            l = l+1;
        end
    else
        solvingOrderSubProblem_i = solvingOrderPreviousSubProblems(solvingOrderPreviousSubProblems(:,5)==solvingOrderPreviousSubProblems(l,5), :);
        
        [rowSubProblem, ~] = find(solvingOrderPreviousSubProblems(l:end,5)>solvingOrderSubProblem_i(1,5),1,'first');
        if isempty(rowSubProblem)
            l = size(solvingOrder,1)+1;
        else
            l = l+rowSubProblem-1;
        end

        [xpInit, yInit, ~, ~, EqXUZ] = sim.solveSubProblemFixZ(solvingOrderSubProblem_i, EqXUZ, xpInit, yInit, y0, zInit);
    end
end

if ~isempty(ie)
    %% Combination of z for the events sub problem
    %varIdx = solvingOrderSubProblem(:,1);
    varIdx = [];
    for k = 1:length(subProblem)
        varIdx = [varIdx; solvingOrderSubProblem{k}(:,1)];
    end
    zIdx = varIdx(varIdx>sim.sys.n+sim.sys.p);

    nZ = length(zIdx);
    if nZ > 0
        nCombinations = 2^nZ;
        zCombinations = [0;1];
        for m = 2:nZ
            zCombinations = [kron(zCombinations, [1;1]), kron(ones(2^(m-1),1), [0;1])];
        end
    else
        nCombinations = 0;
        zCombinations = [];
    end

    % Sort combinations by the number of different z values to zguess
    zPrevious = zInit(zIdx-sim.sys.n-sim.sys.p);
    [~,sortCombinations] = sort(sum(zPrevious~=zCombinations,2),'descend');
    zCombinations = logical(zCombinations(sortCombinations,:));
    
    %% Go through all combinations until a correct solution is found (equalities, inequalites and negative total derivative of event is met)
    EqXUZinit = EqXUZ;
    
    k = 1;
    foundSolution = 0; 
    while (k <= nCombinations) && (foundSolution == 0)
        %% Solve events sub-problem for certain combination of z
        EqXUZ = EqXUZinit;
    
        xpI = xpInit;
        yI = yInit;
        zI = zInit;
        zI(zIdx-sim.sys.n-sim.sys.p) = zCombinations(k,:);
    
        stopMarker = 0;
        if solvingOrderSubProblem{1}(1,3) == 1
            n = solvingOrderSubProblem{1}(1,1);
            v = zCombinations(k,1);
    
            EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
            
        else
            [xpI, yI, G, ~, EqXUZ] = sim.solveSubProblemFixZ(solvingOrderSubProblem{1}, EqXUZ, xpI, yI, y0, zI);
        end
        %% Calculate the solution for the equations of the following sub-Problems from the same sub-set (xp, y for the certain combinations of z)
        l = 1;
        while (l <= size(solvingOrderFollowingSubProblems,1)) %&& (stopMarker==0)
            if solvingOrderFollowingSubProblems(l, 3) == 1
                if solvingOrderFollowingSubProblems(l,1) <= sim.sys.n+sim.sys.p
                    n = solvingOrderFollowingSubProblems(l,1);
                    e = solvingOrderFollowingSubProblems(l,2);
        
                    v = explicitSolution(Hphi(e,:), EqXUZ, HFc(n,:), HFt(n,:), HFf(n,:));
        
                    EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
                    l = l+1;
        
                    if n<=sim.sys.n
                        xpI(1,n) = v;
                    else
                        yI(1,n - sim.sys.n) = v;
                    end
                else
                    n = solvingOrderFollowingSubProblems(l,1);
                    e = solvingOrderFollowingSubProblems(l,2);
        
                    HphiXU = Hphi(e,:);
                    zguess = z0(1, n - sim.sys.n - sim.sys.p);
                    EqXUZ0 = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(zguess-1) - HFf(n,:)'*zguess) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*zguess);
                    F0 = HphiXU * EqXUZ0;
    
                    if sign(F0) <= 0
                        v = zguess;
                    else
                        v = ~zguess;
                    end
                    EqXUZ = EqXUZ .* (1 + HFt(n,:)'.*HFf(n,:)' + HFt(n,:)'*(v-1) - HFf(n,:)'*v) .* (1 - abs(HFc(n,:)') + HFc(n,:)'*v);
                    zI(1,n - sim.sys.n - sim.sys.p) = v;
                    l = l+1;
                end
            else
                solvingOrderSubProblem_i = solvingOrderFollowingSubProblems(solvingOrderFollowingSubProblems(:,5)==solvingOrderFollowingSubProblems(l,5), :);
                [rowSubProblem, ~] = find(solvingOrderFollowingSubProblems(l:end,5)>solvingOrderSubProblem_i(1,5),1,'first');
                if isempty(rowSubProblem)
                    l = size(solvingOrder,1)+1;
                else
                    l = l+rowSubProblem-1;
                end
                
                
                if ismember(solvingOrderSubProblem_i(1,5), subProblem)
                    % sP = find((subProblem == solvingOrderSubProblem_i(1,5)));
                    [~,sP] = ismember(solvingOrderSubProblem_i(1,5), subProblem);
                    [xpI, yI, ~, ~, EqXUZ] = sim.solveSubProblemFixZ(solvingOrderSubProblem{sP}, EqXUZ, xpI, yI, y0, zI);
                elseif solvingOrderSubProblem_i(1,4) == subSet
                    [xpI, yI, zI, ~, ~, EqXUZ] = sim.solveSubProblemVariableZ(solvingOrderSubProblem_i, EqXUZ, xpI, yI, y0, zI);
                else
                    [xpI, yI, ~, ~, EqXUZ] = sim.solveSubProblemFixZ(solvingOrderSubProblem_i, EqXUZ, xpI, yI, y0, zI);
                end
            end
        end
        
        G = HphiIneq * EqXUZ;
        % if k == 1
        %     GInit = G;
        % end
    
        if sum(G>opt.InequalityPositivityTolerance)==0 %1e-5% VERY CRITICAL PARAMETER FOR SIMULATION RESULTS! 5e-9 %(stopMarker == 0) && 
            %% Calculate Derivatives xpp and yp of xp and y        
            J = sim.sys.jacobian(xpI, x0, u0, yI, zI);
            JFxu = [J.equality.state, J.equality.input];
            JFdxy = [J.equality.stateDerivative, J.equality.algebraic];
        
            ddxdy = -JFdxy\(JFxu*[xpI, du]');
            xpp = ddxdy(1:sim.sys.n);
            yp = ddxdy(sim.sys.n+1:end);
        
            %% Check if the total Derivative of the event is negative
            
            % VERY CRITICAL PARAMETER FOR SIMULATION RESULTS! -1e-9
            je = G>=-opt.InequalityZeroTolerance;
            %je = G>=-1e-9;
            %je = ( (G>=-1e-9) & (G/-min(G))>-1e-3 ) | G>=0;

            %je = G>0;
            %[~,idx] = max(G); %isoutlier(G, "ThresholdFactor", 10);%G>-5e-9; %G>-5e-6   
            %je(idx) = true;
            
            %disp(mean(G)+2*std(G))
            %je(ie) = 1;
    
            JG = [J.inequality.stateDerivative, J.inequality.state, J.inequality.algebraic, J.inequality.input];
            DGj = JG * [xpp; xpI'; yp; du'];
            
            if sum(DGj(je)>0)==0
                foundSolution = 1;
                z1 = zI;
                xp1 = xpI;
                y1 = yI;
                dy1 = yp';
            % else
            %     fprintf('!\n')
            end
        end
        k = k+1;
    end
    
    if foundSolution == 0
    warning('continuousSim:NoSolutionForDiscreteEvent','No Solution found for the discrete event!')
    z1 = NaN;
    xp1 = NaN;
    y1 = NaN;
    dy1 = NaN;
    end
else
    J = sim.sys.jacobian(xpInit, x0, u0, yInit, zInit);
    JFxu = [J.equality.state, J.equality.input];
    JFdxy = [J.equality.stateDerivative, J.equality.algebraic];
        
    ddxdy = -JFdxy\(JFxu*[xpInit, du]');
    %xpp = ddxdy(1:sim.sys.n);
    yp = ddxdy(sim.sys.n+1:end);
    if isnan(yp)
        yp(:) = 0;
    end

    z1 = zInit;
    xp1 = xpInit;
    y1 = yInit;
    dy1 = yp';
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

