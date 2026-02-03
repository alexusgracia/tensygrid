function [sim, xsim, ysim, tsim, zsim] = discreteSim(sim, x0, yguess, zguess, ut, u, opt)
% Torben Warnecke - 11/06/2024

[solvingOrder, ~] = sim.determineParallizedSolvingOrder(0); % I is the incidence matrix
%solvingOrderWithoutZ = sim.determineParallizedSolvingOrder(1);

if opt.equationReordering == false
    solvingOrder(:,3) = 0;
    solvingOrder(:,4) = 1;
    solvingOrder(:,5) = 1;
end

sim.solvingOrder = solvingOrder;

%fprintf('Seperation in %d subset(s) with %d subproblem(s), from which %d are explicit solvable.', max(solvingOrder(:,4)), max(solvingOrder(:,5)), sum(solvingOrder(:,3)))


if length(ut) == 2 && isempty(u)
    ut = ut(1):sim.sys.ts:ut(2);
elseif ut(2)-ut(1)  ~= sim.sys.ts
    %     error('Step size of ut does not equal step size ts of the model!')
end

if ~isempty(u) && (size(u,1) ~= length(ut))
    error('Length of input vectors must be similar as the length of the time vector!')
end

xsim = nan([length(ut), sim.sys.n]);
ysim = nan([length(ut),sim.sys.p]);
tsim = ut;
zsim = nan([length(ut),sim.sys.q]);
F = nan([length(ut),sim.sys.nEq]);
G = nan([length(ut),sim.sys.nIneq]);

%try
if isempty(yguess)
    yguess0 = zeros([1, sim.sys.p]);
else 
    yguess0 = yguess;
end
if isempty(zguess)
    zguess = ones([1, sim.sys.q]); %zeros
end
xpguess = x0;

if isempty(u)
    [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, [], solvingOrder, yguess0, xpguess, zguess);
else
    if sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
        [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, u(1,:), sim.solvingOrder, yguess0, xpguess, zguess);
    elseif sim.sys.discreteDifferentiation == "backward"
        [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, u(2,:), sim.solvingOrder, yguess0, xpguess, zguess);
    end
    %[xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, u(1,:), solvingOrder, yguess, xpguess, zguess);
end
if any(isnan(Fi)) || any(isnan(Gi)) || any(abs(Fi)>eps) || any(Gi>eps)
    yguess0 = ones([1, sim.sys.p]);
    zguess = zeros([1, sim.sys.q]);
    xpguess = x0;

    if isempty(u)
        [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, [], solvingOrder, yguess0, xpguess, zguess);
    else
        if sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
            [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, u(1,:), sim.solvingOrder, yguess0, xpguess, zguess);
        elseif sim.sys.discreteDifferentiation == "backward"
            [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, u(2,:), sim.solvingOrder, yguess0, xpguess, zguess);
        end
        %[xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, u(1,:), solvingOrder, yguess, xpguess, zguess);
    end
    if (norm(FiNew) <= norm(Fi)) && (norm(GiNew(GiNew>0)) <= norm(Gi(Gi>0)))
        Fi = FiNew;
        Gi = GiNew;
        xp = xpNew;
        y = yNew;
        z = zNew;
    end
end

xsim(1,:) = x0;
xsim(2,:) = xp;
if sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
    ysim(1,:) = y;
    zsim(1,:) = z;
elseif sim.sys.discreteDifferentiation == "backward"
    ysim(2,:) = y;
    zsim(2,:) = z;
end

F(1,:) = Fi;
G(1,:) = Gi;

% if any(Gi>0)
%     disp('!')
% end

k=1;
while ((sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none") && k < length(ut)) ...
        || (sim.sys.discreteDifferentiation == "backward" && k+1 < length(ut))

    if ~any(isnan(xsim(k,:))) && ...
            ~( sim.sys.discreteDifferentiation ~= "backward" && any(isnan(ysim(k,:))) || sim.sys.discreteDifferentiation == "backward" && any(isnan(ysim(k+1,:)))) ...
            && ~( sim.sys.discreteDifferentiation ~= "backward" && any(isnan(zsim(k,:))) || sim.sys.discreteDifferentiation == "backward" && any(isnan(zsim(k+1,:))))

        if sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
            xsim(k+1,:) = xp;
        end

        k = k+1;

        xpguess = xp;
        yguess = y;
        x0 = xp;

        % try
        
        if isempty(u)
            [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, [], sim.solvingOrder, yguess, xpguess, z);
        elseif sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
            [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, u(k,:), sim.solvingOrder, yguess, xpguess, z);
        elseif sim.sys.discreteDifferentiation == "backward"
            [xp, y, z, Fi, Gi] = sim.solveStepBruteForce(x0, u(k+1,:), sim.solvingOrder, yguess, xpguess, z);
        end
        if any(isnan(Fi)) || any(isnan(Gi)) || any(abs(Fi)>eps^(1/4)) || any(Gi>eps)
            yguessNew = yguess0;
            if isempty(u)
                [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, [], sim.solvingOrder, yguessNew, xpguess, z);
            elseif sim.sys.discreteDifferentiation == "forward" || sim.sys.discreteDifferentiation == "none"
                [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, u(k,:), sim.solvingOrder, yguessNew, xpguess, z);
            elseif sim.sys.discreteDifferentiation == "backward"
                [xpNew, yNew, zNew, FiNew, GiNew] = sim.solveStepBruteForce(x0, u(k+1,:), sim.solvingOrder, yguessNew, xpguess, z);
            end
            if (norm(FiNew) <= norm(Fi)) && (norm(GiNew(GiNew>0)) <= norm(Gi(Gi>0)))
                Fi = FiNew;
                Gi = GiNew;
                xp = xpNew;
                y = yNew;
                z = zNew;
            end
        end
        % catch
        %     disp('!')
        % end
        
        if sim.sys.discreteDifferentiation == "backward"
            xsim(k+1,:) = xp;
            ysim(k+1,:) = y;
            zsim(k+1,:) = z;
        else
            ysim(k,:) = y;
            zsim(k,:) = z;
        end

        F(k,:) = Fi;
        G(k,:) = Gi;
    else
        xsim(k+1:end,:) = nan;
        ysim(k+1:end,:) = nan;
        zsim(k+1:end,:) = nan;
        %tsim = ut(1:k-1);
        warning('discreteSim:SimulationStopped', sprintf("Simulation stopped at t=%f, since no feasible solution could be found for the next time step.", tsim(k)))
        k = length(ut);
    end
end

sim.x = xsim;
sim.y = ysim;
sim.tsim = tsim;
sim.z = zsim;
sim.zt = tsim;

sim.F = F;
sim.G = G;

if max(abs(sim.F),[],'all')>eps^(1/4)
    warning('discreteSim:Inaccurate',"Results may be inaccurate. Absolute function value is bigger then eps^(1/4).\n"...
        + sprintf(['Maximum absolute function value is %0.4g\n'],max(abs(F),[],'all'))...
        + sprintf('Mean absolute function value is %0.4g\n',mean(abs(F),"all")))
end

end
