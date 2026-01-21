classdef dmsim
    % DMSIM Create dmsim object and simulate time response of
    % dynamic (discrete- or continuous-time) implicit MTI models
    % to arbitrary inputs.
    %
    %   input parameter:
    %   - sys: The implicit MTI model stored as a dmss(...) object.
    %   - x0: Initial state values for simulating a state-space
    %   model, specified as a vector having one entry for each
    %   state in sim.sys.
    %   - ut: Time samples at which to compute the response.
    %   In discrete-time simulations t is specified as a vector of
    %   the form T0:dT:Tf. The time step dT must equal the
    %   sample time of the discrete-time simMTI.sys.
    %   In continuous-time simulations the timesamples correspond
    %   to the time where an input value of u is located.
    %   - u: input signals for simulation, u is an array with as
    %   many rows as there are time samples (length(t)) and as many
    %   columns as there are inputs to the system sim.sys.
    %   - opt: solver options, e.g. created with
    %   options = odeset(Name,Value,...).
    %   Options are only available for continuous-time simulations.
    %
    %   properties:
    %   - sys: The implicit MTI model stored as a dmss(...) object.
    %   - ut: time signal with respect to the input signals
    %   - u: input signals
    %   - x: simulated state signals
    %   - tsim: time signal with respect to the simulation results (state
    %   and algebraic signals)
    %   - y: simulated algebraic signals (e.g. output signals)
    %
    % For detailed documentation see <a href="matlab:open((which('dmsimDoc.html')))">here</a>

    % Torben Warnecke - 11/06/2024


    properties
        % MODEL
        sys dmss
        
        % SIMULATION INPUT
        ut
        u
        
        % SIMULATION OUTPUT
        x
        tsim
        y
        
        % BINARY SIMULATION OUTPUT
        z
        zt
        
        % FUNCTION VALUES
        F
        G
    end
    properties (Hidden)
        solvingOrder
        doingParameterID logical = 0

        
    end

    methods
        function sim = dmsim(sys, x0, y0, z0, ut, u, opt)%, maxOdeTime, inequalityTolerance)
            % DMSIM Create dmsim object and simulate time response of
            % dynamic (discrete- or continuous-time) implicit MTI models
            % to arbitrary inputs.
            %
            %   input parameter:
            %   - sys: The implicit MTI model stored as a dmss(...) object.
            %   - x0: Initial state values for simulating a state-space
            %   model, specified as a vector having one entry for each
            %   state in sim.sys.
            %   - y0: Initial guess for the initial algebraic values. If y0
            %   is empty (y0 = []), the initial guess is zero. If it fails
            %   the function will try with ones.
            %   -z0: Initial guess for the initial boolean values. If z0
            %   is empty (z0 = []), the initial guess is zero.
            %   - ut: Time samples at which to compute the response.
            %   In discrete-time simulations t is specified as a vector of
            %   the form T0:dT:Tf. The time step dT must equal the
            %   sample time of the discrete-time simMTI.sys.
            %   In continuous-time simulations the timesamples correspond
            %   to the time where an input value of u is located.
            %   - u: input signals for simulation, u is an array with as
            %   many rows as there are time samples (length(t)) and as many
            %   columns as there are inputs to the system sim.sys.
            %   - opt: solver options, e.g. created with
            %   options = odeset(Name,Value,...).
            %   Options are only available for continuous-time simulations.
            %   Additionally the options: InequalityZeroTolerance,
            %   InequalityPositivityTolerance and MaxTime have been added,
            %   which can be changed/added by e.g.
            %      opt = odeset('RelTol', 1e-9, 'AbsTol', 1e-3);
            %      opt.InequalityZeroTolerance = 1e-6;
            %      opt.InequalityPositivityTolerance = 1e-3;
            %      opt.MaxTime = inf;
            %
            %   properties:
            %   - sys: The implicit MTI model stored as a dmss(...) object.
            %   - ut: time signal with respect to the input signals
            %   - u: input signals
            %   - x: simulated state signals
            %   - tsim: time signal with respect to the simulation results (state
            %   and algebraic signals)
            %   - y: simulated algebraic signals (e.g. output signals)
            %
            % For detailed documentation see <a href="matlab:open((which('dmsimDoc.html')))">here</a>

            % Torben Warnecke - 11/06/2024

            sim.sys = sys;
            sim.sys.H = sim.sys.H.convert2doubles();
            sim.sys = sim.sys.normalizeParameters(1e0);

            if nargin > 1
                if length(x0) ~= sim.sys.n
                    error('Number of initial states does not equal the number of states from the model!')
                elseif sim.sys.n>0 && size(x0, 1) ~= 1
                    x0 = x0';
                    if size(x0, 1) ~= 1
                        error('Initial states x0 must be parsed as vector.')
                    end
                end
                if nargin > 5
                    if size(u,2) ~= sim.sys.m
                        fprintf('Number of inputs: %d\n', size(u,2))
                        fprintf('Number of inputs necessary: %d\n', sim.sys.m)
                        error('Number of inputs does not equal the number of inputs from the model! Try transposing the input-vector/matrice.')
                    end
                end
                if nargin < 2
                    error('Not enough arguments. Initial state values are needed!')
                else
                    x0 = reshape(x0, [], length(x0));
                end
                
                if nargin < 3
                    y0 = zeros([1, sim.sys.p]);
                elseif (length(y0) ~= sim.sys.p) && (~isempty(y0))
                    error(['Number of initial algebraic variable guesses does not equal the number of algebraic variables from the model!' ...
                        'Either must be equal to sys.p or empty (y0 = [])!'])
                end
                if nargin < 4
                    z0 = zeros([1, sim.sys.p]);
                elseif (length(z0) ~= sim.sys.q) && (~isempty(z0))
                    error(['Number of initial boolean variable guesses does not equal the number of boolean variables from the model!' ...
                        'Either must be equal to sys.q or empty (z0 = [])!'])
                end

                if nargin < 5
                    if sim.sys.ts == 0
                        ut = [0, 10];
                    else
                        ut = 0:sim.sys.ts:10;
                    end
                end
                sim.ut = ut;

                % if nargin <9
                %     inequalityTolerance = 1e-9;
                % end

%                 if license('test', 'Distrib_Computing_Toolbox')
%                     p = gcp('nocreate');
%                     if ~isempty(p)
%                         delete(gcp)
%                     end
%                     parpool('local')
%                 end

                if sim.sys.ts == 0
                    if nargin < 7
                        %opt = odeset('RelTol', 10.0^(-3),'AbsTol',10.0^(-7));
                        opt = odeset();
                    end
                    % if nargin < 8
                    %     maxOdeTime = inf;
                    % elseif isempty(maxOdeTime)
                    %     maxOdeTime = inf;
                    % end
                    if ~isfield(opt, 'InequalityZeroTolerance')
                        opt.InequalityZeroTolerance = 1e-6; %1e-9
                    end
                    if ~isfield(opt, 'InequalityPositivityTolerance')
                        opt.InequalityPositivityTolerance = 1e-3; %1e-5
                    end
                    if ~isfield(opt, 'MaxTime')
                        opt.MaxTime = inf;
                    end


                    if nargin > 5
                        sim.u = u;
                        [sim, ~, ~, ~, ~, ~] = sim.continuousSimTotalDerivative(x0, y0, z0, ut, u, opt); %, maxOdeTime, inequalityTolerance);
                    else
                        sim.u = [];
                        u = [];
                        [sim, ~, ~, ~, ~, ~] = sim.continuousSimTotalDerivative(x0, y0, z0, ut, u, opt); %, maxOdeTime, inequalityTolerance);
                    end
                else
                    if nargin < 7
                        opt.equationReordering = true;
                    elseif ~isfield(opt, "equationReordering")
                        opt.equationReordering = true;
                    elseif isfield(opt, "equationReordering")
                        if opt.equationReordering ~= true && opt.equationReordering ~= false 
                            error("opt.equationReordering can be either true or false. Default is true.")
                        end
                    end
                    if nargin > 5
                        sim.u = u;
                        [sim, ~, ~, ~, ~] = sim.discreteSim(x0, y0, z0, ut, u, opt);
                    else
                        sim.u = [];
                        u = [];
                        [sim, ~, ~, ~, ~] = sim.discreteSim(x0, y0, z0, ut, u, opt);
                    end
                end
            end
        end
        function sim = set.sys(sim, newSys)
            sim.sys = newSys;
            if sim.sys.n+sim.sys.p+sim.sys.q ~= sim.sys.nIneq + sim.sys.nEq
                error('Not a square system! Number of equation and inequality constraints are not equal to the number of state, algebraic and boolean variables.')
            end
        end
    end
end