classdef simMTI < handle
    properties
        sys mti
    end

    properties (Access = public)
        x0 % simulation initial states
        tsim % simulation time vector
        xsim % simulated state trajectory
        ysim % simulated output trajectory
        usim % simulation input vector
    end

    methods
        function sim = simMTI(msys,varargin)
            sim.sys = msys;

            if(nargin >=2)
                error('Too many input arguments');
            end
        end
        
        function plot(thisSim)
            %PLOTSIM Plot simulation results
            % The simulation results are summarized in one figure. The
            % inputs as well as the state and output trajectories are
            % plotted over time.
            % INPUTS:
            % ---
            % OUTPUTS:
            % ---
            figure
            subplot(2,2,1)
            if thisSim.sys.m == 0
                text(0.05,0.5,'No inputs','FontSize',14)
            else
                plot(thisSim.tsim,thisSim.usim)
                legend(thisSim.createLabel('u',thisSim.sys.m))
            end
            title('Inputs')
            subplot(2,2,2)
            plot(thisSim.tsim,thisSim.xsim)
            legend(thisSim.createLabel('x',thisSim.sys.n))
            title('States')
            subplot(2,2,[3 4])
            plot(thisSim.tsim,thisSim.ysim)
            legend(thisSim.createLabel('y',thisSim.sys.p))
            title('Outputs')
        end

        function [y,tOut,x] = simulate(sim, u, t, x0)
            arguments
                sim
                u
                t
                x0
            end
            % SIMULATE - simulate time response of dynamic (discrete- or continuous-time) explicit MTI models to arbitrary inputs. 
            %
            %   input parameter:
            %   - sim: contains the MTI simulation object simMTI(...) with
            %   the according system stored as a mss(...) object
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
            %   sim = simMTI(msys);
            %
            %   x0 = [-0.2 0.3];
            %   t = 0:0.05:8;
            %   u = zeros(length(1),1)';
            %   u(t>=2) = 1;
            %
            %   [y,tOut,x] = sim.simulate(u, t, x0,includeStateTrajectory);
            %
            
            % Enrico Uhlenberg , Torben Warnecke, Carlos Cateriano
            % Yáñez, Leona Schnelle - 12/06/2024


            % Check if system is time-discrete or time-continuous
            if sim.sys.ts > 0
                % Discrete-time Simulation
                if(nargout==1)
                    [y]= msimDisc(sim, u, t, x0);
                elseif (nargout ==2)
                    [y,tOut]= msimDisc(sim, u, t, x0);
                elseif (nargout ==3)
                    [y,tOut,x]= msimDisc(sim, u, t, x0);
                end
            elseif sim.sys.ts == 0
                % Continuous-time Simulation
                method = 'linear'; % we may implement different integration methods for the inputs, if possible
                if(nargout==1)
                    [y]= msimCont(sim, u, t, x0, method);
                elseif (nargout ==2)
                    [y,tOut]= msimCont(sim, u, t, x0, method);
                elseif (nargout ==3 || nargout ==0)
                    [y,tOut,x]= msimCont(sim, u, t, x0, method);
                end
            else
                error('Unallowed step size ts, the step size needs to be zero (continuous-time) or bigger then zero (discrete-time).')
            end
            
            if nargout == 0
                sim.tsim = tOut;
            end
            
            if nargout > 2
                sim.xsim = x;
            end
            
            sim.ysim = y; % TODO Is the output vector one elemnt shorter than the statevector? The old msim seems to indicate that
            sim.usim = u;
        end
    end

    methods (Access = protected)

        function labels = createLabel(~,plot_var,var_number)
            %CREATELABEL Create legend entries for plots
            % Legend entries for the plots of the simulation results are
            % build for inputs, states and outputs
            % INPUTS:
            % - plot_var: name of the variable ('x', 'u' or 'y')
            % - var_numbers: number of legend entries for this variable
            % OUTPUTS:
            % - labels: Cell array with legend entries
            labels = cell(1,var_number);
            for ind1 = 1:var_number
                labels{ind1} = [plot_var num2str(ind1)];
            end
        end
    end
end