function [sim, xsim, ysim, tsim, zsim, zt] = continuousSimTotalDerivative(sim, x0, yguess, zguess, ut, u, opt)%, maxOdeTime, inequalityTolerance)
% Torben Warnecke - 11/06/2024

[solvingOrder, ~] = sim.determineParallizedSolvingOrder(0); % I is the incidence matrix
%solvingOrderWithoutZ = sim.determineParallizedSolvingOrder(1);

% solvingOrder(:,3) = 0;
% solvingOrder:,4) = 1;
% solvingOrder(:,5) = 1;

sim.solvingOrder = solvingOrder;

%fprintf('Seperation in %d subset(s) with %d subproblem(s), from which %d are explicit solvable.', max(solvingOrder(:,4)), max(solvingOrder(:,5)), sum(solvingOrder(:,3)))

tsim = [];
xsim = [];
ysim = [];
zsim = [];
zt = [ut(1)];

t = ut(1);

te = ut(1);

% try
if isempty(yguess)
    yguess = zeros([1, sim.sys.p]);
end
if isempty(zguess)
    zguess = zeros([1, sim.sys.q]);
end

    ie = solvingOrder(find(solvingOrder(:,2)>sim.sys.nEq, 1), 2) - sim.sys.nEq;
    [xp0, y0, z, dy0, Fi, Gi] = sim.eventTotalDerivative(ie, ut(1), x0, yguess, zguess, u, ut, solvingOrder, opt);
    % if any(isnan(xp0))
    %     error('No solution found!')
    % end
% catch


if any(isnan(xp0)) || any(isnan(dy0)) || any(isnan(Fi)) || any(isnan(Gi)) || any(abs(Fi)>opt.InequalityZeroTolerance) || any(Gi>opt.InequalityPositivityTolerance)
    yguess = ones([1, sim.sys.p]);
    zguess = zeros([1, sim.sys.q]);

    ie = solvingOrder(find(solvingOrder(:,2)>sim.sys.nEq, 1), 2) - sim.sys.nEq;
    [xp0New, y0New, zNew, dy0New, FiNew, GiNew] = sim.eventTotalDerivative(ie, ut(1), x0, yguess, zguess, u, ut, solvingOrder, opt);
    if any(isnan(xp0New)) || any(isnan(dy0New))
        warning('continuousSim:NoInitialBinaryVarFound', sprintf("Simulation stopped at ut(1), since no feasible initial solution could be found."))
        xsim = nan([2, sim.sys.n]);
        ysim = nan([2, sim.sys.p]);
        zsim = nan([2, sim.sys.q]);
        t = [ut(1); ut(end)];
        tsim = t;
        zt = [ut(1); ut(end)];
        te = [];
    elseif (norm(FiNew) <= norm(Fi)) && (norm(GiNew(GiNew>0)) <= norm(Gi(Gi>0)))
        xp0 = xp0New;
        y0 = y0New;
        z = zNew;
        dy0 = dy0New;
        Fi = FiNew;
        Gi = GiNew;
    elseif any(isnan(xp0)) || any(isnan(dy0)) || any(isnan(Fi)) || any(isnan(Gi))
        xp0 = xp0New;
        y0 = y0New;
        z = zNew;
        dy0 = dy0New;
        Fi = FiNew;
        Gi = GiNew;
    end
end

if ~any(isnan([xsim, ysim, zsim]))
    zsim = z;
    var0 = [x0, y0];
    dvar0 = [xp0, dy0];
end
tOde = tic();

while ~isempty(te) && t(end)~=ut(end)
    % ODE Set up
    if license('test', 'Symbolic_Toolbox') && 1 == 0
        % handle = symbolicFunctionHandleGenerator(sim.sys, ut, u);
        % 
        % udae = @(t,k)interpolatedInput(sim.sys,t,ut,u(:,k));
        % binaries = reshape(z, 1, sim.sys.q);
        % 
        % if sim.sys.q>0
        %     handle = @(t,x,xp) iMTI_as_DAE_Temp(t,x,xp,binaries,arrayfun(udae,t*ones(1,sim.sys.m), 1:sim.sys.m));
        % else
        %     handle = @(t,x,xp) iMTI_as_DAE_Temp(t,x,xp,arrayfun(udae,t*ones(1,sim.sys.m), 1:sim.sys.m));
        % end
        % 
        % if sim.sys.m > 0
        %     eventoption = odeset('Events',@(t,x,xp) inequalityEvent(sim.sys,t,x,xp,z,ut,u,maxOdeTime));
        % else
        %     eventoption = odeset('Events',@(t,x,xp) inequalityEvent(sim.sys,t,x,xp,z,ut,[],maxOdeTime));
        % end
    else
        if sim.sys.m > 0
            handle = @(t,x,xp) odeHandle(sim.sys, t, x, xp, z, ut, u);
            eventoption = odeset('Events',@(t,x,xp) inequalityEvent(sim.sys,t,x,xp,z,ut,u, opt.MaxTime));
            jacobioption = odeset('Jacobian', @(t,x,xp) jacobiHandle(sim.sys,t,x,xp,z,ut,u));
        else
            handle = @(t,x,xp) odeHandle(sim.sys, t, x, xp, z, ut,[]);
            eventoption = odeset('Events',@(t,x,xp) inequalityEvent(sim.sys,t,x,xp,z,ut,[], opt.MaxTime));
            jacobioption = odeset('Jacobian', @(t,x,xp) jacobiHandle(sim.sys,t,x,xp,z,ut,[]));
        end
        opt.Jacobian = jacobioption.Jacobian;
    end
    
    opt.Events = eventoption.Events;
    

    % simulate until inequality breach
    tspan = [te(end) ut(end)];

    [t,x,te,xe,ie] = ode15i(handle,tspan,var0',dvar0', opt);
    tsim = [tsim; t];
    xsim = [xsim; x(:, 1:sim.sys.n)];
    ysim = [ysim; x(:, sim.sys.n+1:sim.sys.n+sim.sys.p)];
    zsim = [zsim; z];

    if ~isempty(ie) && any(ie > sim.sys.nIneq)
        % Calculation Time of ode solver exceeded maxOdeTime given by user.
        %zsim = [zsim; z];
        te = [];
        warning('continuousSim:CalculationTimeExceeded', 'Calculation Time of the ode solver exceeded maxOdeTime given by user! Simulation stopped. Results can be inaccurate.')

    elseif ~isempty(te) && (te(end) > tspan(1)) && ~ismember(te(end), zt) && (t(end)==te(end))
        % inequality boundary is hit!
        zt = [zt; te(end)];

        t0 = tsim(end);
        x0 = xsim(end,:);
        y0 = ysim(end,:);
        z0 = z;

        dx0 = (x0-xsim(end-1,:))/(t0-tsim(end-1));
        for l = 1:size(u,2)
           u0(1,l) = interp1q(ut, u(:,l), t0);
        end

        [dx0, y0, z0, Fe, Ge] = sim.solveStepBruteForce(x0, u0, solvingOrder, y0, dx0, z0);
        %je = find((Ge/-min(Ge))>-1e-3 | Ge>=0);
        %je = find(Ge>=-1e-9);
        je = find(Ge>=-opt.InequalityZeroTolerance);
        %find(((Ge>=-1e-9) & (Ge/-min(Ge))>-1e-3 ) | Ge>=0);
        je = unique([je;ie]);

        %[xp1, ~, z1, dy1] = sim.eventTotalDerivative(ie(end), t0, x0, y0, z0, u, ut, solvingOrder);
        % [xp1, ~, z1, dy1] = sim.eventTotalDerivative(je, t0, x0, y0, z0, u, ut, solvingOrder, opt);
        [xp1, y1, z1, dy1] = sim.eventTotalDerivative(je, t0, x0, y0, z0, u, ut, solvingOrder, opt);
        %disp(xp1)

        if any(isnan(z1)) || any(isnan(xp1))
            warning('continuousSim:NoSolutionForDiscreteEvent','No Solution found for the discrete event!')
            te = [];
            zsim = [zsim; z];
        else
            z = z1;
            zsim = [zsim; z];
            zt = [zt; te(end)];
            %var0 = [xe(end,:)];
            var0 = [x0, y1];
            dvar0 = [xp1, dy1];
        end

    elseif ~isempty(te) && (te(end) == tspan(1)) %&& ismember(te(end), zt) && (t(end)==te(end))
        warning('continuousSim:ZenoBehaviuour','Zeno behaviour detected. Simulation stopped. Implementing a low-pass filter might be a solution for simulating the given system.')
        te = [];
        %zsim = [zsim; z];
    end
end

zt = [zt; tsim(end)];

    function F = odeHandle(sys,t,x,xp,z,ut,u)
        if ~isempty(u)
            uode = zeros(sys.m,1);
            for k = 1:size(u,2)
                %uode(k,1) = interp1(ut, u(:,k), t, method);
                uode(k,1) = interp1q(ut, u(:,k), t);
            end
        else
            uode = [];
        end
        
        dxode = xp(1:sys.n);
        xode = x(1:sys.n);
        yode = x(sys.n+1:sys.n+sys.p);

        [Eq, ~, ~] = sys.createEquations(dxode, xode, yode, uode, z);

        F = Eq;
    end

    function [value,isterminal,direction] = inequalityEvent(sys,t,x,xp,z,ut,u, maxOdeTime)
        if ~isempty(u)
            uode = zeros(sys.m,1);
            for k = 1:size(u,2)
                %uode(k,1) = interp1(ut, u(:,k), t, method);
                uode(k,1) = interp1q(ut, u(:,k), t);
            end
        else
            uode = [];
        end

        dxode = xp(1:sys.n);
        xode = x(1:sys.n);
        yode = x(sys.n+1:sys.n+sys.p);

        [~, Ineq, ~] = sys.createEquations(dxode, xode, yode, uode, z);

        tOdeCurrent = toc(tOde);
        value = [Ineq; tOdeCurrent-maxOdeTime]; %- eps^0.5;
        isterminal = ones([sys.nIneq+1,1]);
        direction = ones([sys.nIneq+1,1]);
        %direction = zeros([sys.nIneq,1]);
    end

    function [dfdx, dfdp] = jacobiHandle(sys,t,x,xp,z,ut,u)
        if ~isempty(u)
            uode = zeros(sys.m,1);
            for k = 1:size(u,2)
                %uode(k,1) = interp1(ut, u(:,k), t, method);
                uode(k,1) = interp1q(ut, u(:,k), t);
            end
        else
            uode = [];
        end

        dxode = xp(1:sys.n);
        xode = x(1:sys.n);
        yode = x(sys.n+1:sys.n+sys.p);

        J = sys.jacobian(dxode, xode, uode, yode, z);

        dfdx = [J.equality.state, J.equality.algebraic];
        dfdp = [J.equality.stateDerivative, zeros(sys.n +sys.p, sys.p)];
    end

    % function handle = symbolicFunctionHandleGenerator(sys, ut, u)
    %     t = sym("t");
    %     vars = [];
    %     xDAE = [];
    %     y = [];
    %     uDAE = [];
    %     for k = 1:sys.n
    %         xDAE = [xDAE; symfun(sprintf('x%u(t)',k),t)];
    %     end
    %     for k = 1:sys.p
    %         y = [y; symfun(sprintf('y%u(t)',k),t)];
    %     end
    %     uDAE = [];
    %     for k = 1:sys.m
    %         uDAE = [uDAE, str2sym(sprintf('u%u(t)',k))];
    %     end
    %     zAsDaeParameter = sym('z', [1,sys.q]);
    % 
    %     vars = [xDAE; y];
    %     eqs = sys.createEquations(diff(xDAE, t), xDAE, y, uDAE, zAsDaeParameter);
    % 
    %     if sys.q >0
    %         handle = daeFunction(eqs, formula(vars),zAsDaeParameter, uDAE, 'File','iMTI_as_DAE_Temp');
    %     else
    %         handle = daeFunction(eqs, formula(vars), uDAE, 'File','iMTI_as_DAE_Temp');
    %     end
    % end
    
    function uode = interpolatedInput(sys, t, ut, u)
        if ~isempty(u)
            uode = interp1q(ut,u,t);
        else
            uode = [];
        end
    end


sim.x = xsim;
sim.y = ysim;
sim.tsim = tsim;
sim.z = zsim;
sim.zt = zt;
% delete iMTI_as_DAE_Temp.m
end

