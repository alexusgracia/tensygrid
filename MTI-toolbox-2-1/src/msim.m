function [y,tOut,x] = msim(msys,u,t,x0, interpolationMethod)
%MSIM Create simMTI object and simulate time response of dynamic (discrete- or continuous-time) explicit MTI models to arbitrary inputs.
%
%   input parameter:
%   - msys: Either an explicit MTI model stored as a mss(...) object
%   or an implicit MTI model stored as dmss(...) object.
%   - u: input signals for simulation, u is an array with as
%   many rows as there are time samples (length(t)) and as many
%   columns as there are inputs to the system simMTI.sys.
%   - t: Time samples at which to compute the response.
%   In discrete-time simulations t is specified as a vector of
%   the form T0:dT:Tf. The time step dT must equal the
%   sample time of the discrete-time simMTI.sys.
%   In continuous-time simulations the timesamples correspond
%   to the time where an input value of u is located.
%   For discrete-time simMTI.sys,
%   - x0: Initial state values for simulating a state-space
%   model, specified as a vector having one entry for each
%   state in simMTI.sys.
%   - interpolationMethod: specifies the interpolation method used for the
%   simulated output and state signals of continuous time simulation when 
%   t is a sampled time array.
%   For the interpolation Method "none" (default) the timesteps of the ode solver 
%   will be returned. For "zoh" assumes the signals are
%   piecewise constant. For "foh" assumes the signals are piecewise linear.
%
%   output parameter:
%   - y: Simulated response data (outputs), returned as an
%   array. y is an array with as many rows as there are time
%   samples (length(t)) and as many columns as there are
%   outputs in simMTI.sys.
%   - tOut: Time vector used for simulation, returned as a
%   column vector. In discrete-time tOut is equal to the input
%   parameter t. In continuous time tOut is specified by the
%   ODE solver.
%   - x: State trajectories, returned as an array.
%   x is an array with as many rows as there are time samples
%   (length(t)) and as many columns as there are states in
%   simMTI.sys.
%
%   Example: Simulating a LTI model as a eMTI model
%   A = [-1.5 -3; 3 -1];
%   B = [1.3; 0];
%   C = [1.15 2.3];
%   D = [0];
%
%   F = CPN1(diag(ones([3,1])), [A, B]);
%   G = CPN1(diag(ones([3,1])), [C, D]);
%   msys = mss(F,G);
%
%   x0 = [-0.2 0.3];
%   t = 0:0.05:8;
%   u = zeros(length(1),1)';
%   u(t>=2) = 1;
%
%   [y, tOut, x] = msim(msys, u', t, x0);
%
%   For detailed documentation see <a href="matlab:open((which('msimDoc.html')))">here</a>
% Enrico Uhlenberg, Torben Warnecke, Carlos Cateriano Yáñez - 12/06/2024
arguments
    msys
    u
    t (:,1)
    x0
    interpolationMethod = "none"
end

if ~(interpolationMethod == "none" || interpolationMethod == "zoh" || interpolationMethod == "foh")
    error("Wrong input argument for interpolationMethod, must be either none, zoh, foh.")
end

if isa(msys, 'mss')
    osim = simMTI(msys);
    [y,tOut,x]= osim.simulate(u,t,x0);
elseif isa(msys, 'dmss')
    osim = dmsim(msys, x0, t, u);
    y = osim.y;
    tOut = osim.tsim;
    x = osim.x;
else
    error('Wrong input argument for msys, must be either a mss or dmss object.')
end

if ~isequal(t, tOut) && interpolationMethod~="none"
    % If switches occur during simulation, there will be duplicates of
    % time stamps. The second value of these duplicates is the value after
    % the switch. For interp1 no duplicates are allowed.
    [~, idx] = unique(tOut, 'last');
    if interpolationMethod == "zoh"
        y = interp1(tOut(idx), y(idx,:), t', 'previous');
        x = interp1(tOut(idx), x(idx,:), t', 'previous');
        tOut = t;

    elseif interpolationMethod == "foh"
        y = interp1(tOut(idx), y(idx,:), t', 'linear');
        x = interp1(tOut(idx), x(idx,:), t', 'linear');
        tOut = t;

    end
end

end
