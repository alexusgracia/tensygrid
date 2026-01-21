classdef DisaggregatedPipe< dynamicprops
    % Currently does not work with negative mass flow
    % Due to the darcy friction factor. It is too high. 
    % Negative ranges are not yet implemented here !
    % Maybe a i need an interface for the absolute value of m_dot
    % Still open question. How to make absolute value ??
    properties
        di; % pipe inner diameter in m 
        l; % pipe length in meter
        rho=1000; % fluid density in kg / m3
        cp=4190; % fluid specifi heat capacity in J/kgK 
        mu=466*1e-6; % dynamic viscosity
        U=0.1295; % heat transmitance in W/mK --> 1D for pipe
        Tout = 9+273.15; % Outside temperature --> Ground --> constant
        delta_x = 50; % Discretization steps 
        N=2; % Amount of segments --> Minimal of 2
        darcyParamsFile = 'pipeFrictionFactorMTINeu.mat';
        ReLimit = 2300; % Limit of Reynoldsnumber for laminar flow
        usePipeWall = false;
        rhoPipe = 7850; % pipe wall density in kg/m3
        cpPipe = 460; % pipe wall thermal capacity in J/kgK
        pipeThickness = 0.0035; % pipe thickness in m 
        pipeRoughness = 7e-06; % pipe roughness in m
        initiateModel = true;
        node_source = '1';
        node_target = '2';
        pipeName = 'p_1';
        pipeFrictionFactor; % Auxilliary variable for the Darcy Friction factor
        darcy_params; % Approximation parameters for Darcy Friction factor
        pressureModel; % Multilinear pipe pressure loss model
        temperatureModel; % multilinear pipe temperature model
        massModel;
        pipeModel;
    end
    
    methods
        % Constructor of the pipe class
        function obj = DisaggregatedPipe(opts)
            % Enable the possibility of using name-value arguments to set
            % the class property
            arguments
                opts.?DisaggregatedPipe
            end

            assert(~all(isfield(opts,{'N','delta_x'})), ...
                "Cannot insert both segment amount (N)" + ...
                "and segment length (delta_x)")

            % Now replace all the values of the property depending on the
            % given values. If not given, default values are used
            for prop = string(fieldnames(opts))'

                obj.(prop) = opts.(prop);
                switch prop
                    case "N"
                        obj.delta_x = obj.l/obj.N;
                    case "delta_x"
                        obj.N = floor(obj.l/obj.delta_x);
                        if obj.N < 2, obj.N=2; obj.delta_x = obj.l/obj.N; end
                end
            end
            
            if obj.initiateModel
                obj.setDarcyParams();
                % Now generate the pressure model
                obj.TotalPressureModel();
                % Now generate the temperature model
                obj.TotalTemperatureModel();
                % Now generate mass equation
                obj.massEquation();
                % Now connect both models;
                obj.pipeModel = connect(obj.massModel, obj.temperatureModel, obj.pressureModel);
            end
            
        end

        function setDarcyParams(obj)
            
            if ~exist('pipeFrictionFactor','var')
                obj.pipeFrictionFactor = load(obj.darcyParamsFile, 'pipeFrictionFactor');
                obj.pipeFrictionFactor = obj.pipeFrictionFactor.('pipeFrictionFactor');
                obj.darcy_params = obj.pipeFrictionFactor{1,string(obj.pipeRoughness)}{1}{1,string(obj.di)};         
            else
                obj.darcy_params = pipeFrictionFactor{1,string(obj.pipeRoughness)}{1}{1,string(obj.di)};
            end

            %if ~exist('darcyParams','var')
            %    load(obj.darcyParamsFile); 
            %end
            %obj.darcy_params = darcyParams{:,string(obj.di)};
        end

        function booleanLamTurb = booleanLamTurb(obj)
            % Generate a boolean variable (z1) to get the laminar and
            % turbulent region of the flow regime.
            % To connect this MTI model to the other MTI model, z1 is
            % repesented using input (u2) instead.
            % The equation is the following
            % (2*u2-1)*(u1 - md_min) <=0 --> eq(1)
            % u1 = mass flow
            % u2 = "boolean" variable
            
            % Calculate the limit of mass flow for the minimum reynolds
            % number
            mpmin = obj.ReLimit*(pi*obj.di*obj.mu)/4;
            H1 = hyCPN1();
            % The input is in the following repesentation (matrix)
            %       c
            %   u1  0	1	0	1
            %   u2  0	0	1	1
            H1.F.input = [sparse([1,2,1,2],[2,3,4,4],1,2,4)];
            % Directly add it to the inequality
            H1.phi.inequality = [mpmin, -1, -2*mpmin, 2];
            booleanLamTurb = dmss(H1,0);
            booleanLamTurb.inputName = ["mass_flow_abs","Re_bool"];
        end

        function darcyLaminar = darcyLaminar(obj)
            % Calculate the pressure loss on the laminar part
            % 64 = darcyLaminar *  Re
            % Due to this formula, a numerical trick using a small value is
            % required to ensure the calculation works
            % ReFactor * (u1 + eps) * u2 -64 = 0 --> eq (1)
            % u1 = mass flow 
            % u2 = laminar pipe friction factor
            % eps = 1e-6 --> small number for numerical accuracy
            ReFactor = 4/(pi*obj.di*obj.mu);
            H1 = hyCPN1();
            H1.F.input = [sparse([2,2,1],[1,2,2],1,2,3)];
            H1.phi.equality = [ReFactor*(1e-6), ReFactor ,-64];
            darcyLaminar = dmss(H1,0);
            darcyLaminar.inputName = ["mass_flow_abs","lam_friction_factor"];
        end

        function darcyFactor = darcyFactor(obj)
            % Generate an approximation of the darcy friction factor based
            % on the quadratic functions

            % Now put a simple "heaviside" function. A simple change here
            % to ensure x = 0 is put into one
            % myHeav = @(X)double(X>=0);

            % The structure of the problem is the following
            % u1 is the mass incoming mass flow
            % y2 is once again the mass flow 
            % by setting u1 = y2 --> quadratic relations can be made
            % y3 is here the darcy friction factor required 
            % 0 = x1*(x2+x3*u1)*(x4+x5*y2)*(x6+x7*y3) + x8*(x9+x10*u1)*(x11+x12*y2)
            % Now reformulation of to the darcy fiction factor 
            % y3 = ...... --> yields the different x values 
            % in this context called the darcy_params
            % calculated using the lsqcurvefit 
            % options = optimoptions('lsqcurvefit','Algorithm','interior-point','MaxFunctionEvaluations',50e3,'MaxIterations',10e3);

            % from the x values --> structure matrices are generated
            % First we need to norm everything to CPN1
            % example here for the first param for u1
            % x3/(abs(x2)+abs(x3)) * (2*myHeav(x(2)-1)
            % using some number examples 
            % x2 = 0.8620; x3 = 0.7278
            % (0.8620 + 0.7278 * u1) conversion to Norm 1 the following
            % (1-abs(0.4578) + 0.4578 * u1) * (1.5899)

            %  (1.5899) here is then used for the phi matrix !
            % x1 * (abs(x2) + abs(x3) ......
            % s(1,1) = obj.darcy_params(3)/(abs(obj.darcy_params(2))+abs(obj.darcy_params(3)))*(2*myHeav(obj.darcy_params(2))-1);
            % s(2,1) = obj.darcy_params(5)/(abs(obj.darcy_params(4))+abs(obj.darcy_params(5)))*(2*myHeav(obj.darcy_params(4))-1);
            % s(3,1) = obj.darcy_params(7)/(abs(obj.darcy_params(6))+abs(obj.darcy_params(7)))*(2*myHeav(obj.darcy_params(6))-1);
            % s(1,2) = obj.darcy_params(10)/(abs(obj.darcy_params(9))+abs(obj.darcy_params(10)))*(2*myHeav(obj.darcy_params(9))-1);
            % s(2,2) = obj.darcy_params(12)/(abs(obj.darcy_params(11))+abs(obj.darcy_params(12)))*(2*myHeav(obj.darcy_params(11))-1);
            % %s(:,3) = 0;
            % phi(1,1) = obj.darcy_params(1)*(abs(obj.darcy_params(2))+abs(obj.darcy_params(3)))*(abs(obj.darcy_params(4))+abs(obj.darcy_params(5)))*(abs(obj.darcy_params(6))+abs(obj.darcy_params(7)))*(2*myHeav(obj.darcy_params(2))-1)*(2*myHeav(obj.darcy_params(4))-1)*(2*myHeav(obj.darcy_params(6))-1);
            % phi(1,2) = obj.darcy_params(8)*(abs(obj.darcy_params(9))+abs(obj.darcy_params(10)))*(abs(obj.darcy_params(11))+abs(obj.darcy_params(12)))*(2*myHeav(obj.darcy_params(9))-1)*(2*myHeav(obj.darcy_params(11))-1);
            % phi(1,3) = 1e-2;
            
            s = obj.darcy_params.s;
            phi = obj.darcy_params.phi;

            % Now generate the hyCPN1
            H = hyCPN1();
            % To simplify put everything as an input.
            % It is not so problematic as it will be connected to other
            % models on the next steps
            H.F.input = [s(1:3,:)];
            %H.F.algebraic = [s(2:end,:)];
            H.phi.equality = phi;
            darcyFactor = dmss(H,0);
            darcyFactor.inputName = ["mass_flow_abs","mass_flow_square","turb_friction_factor"];
        end

        % function darcyFactor = darcyFactor(obj)
        %     % Generate an approximation of the darcy friction factor based
        %     % on the quadratic functions
        % 
        %     % Now put a simple "heaviside" function. A simple change here
        %     % to ensure x = 0 is put into one
        %     myHeav = @(X)double(X>=0);
        % 
        %     % The structure of the problem is the following
        %     % u1 is the mass incoming mass flow
        %     % y2 is once again the mass flow 
        %     % by setting u1 = y2 --> quadratic relations can be made
        %     % y3 is here the darcy friction factor required 
        %     % 0 = x1*(x2+x3*u1)*(x4+x5*y2)*(x6+x7*y3) + x8*(x9+x10*u1)*(x11+x12*y2)
        %     % Now reformulation of to the darcy fiction factor 
        %     % y3 = ...... --> yields the different x values 
        %     % in this context called the darcy_params
        %     % calculated using the lsqcurvefit 
        %     % options = optimoptions('lsqcurvefit','Algorithm','interior-point','MaxFunctionEvaluations',50e3,'MaxIterations',10e3);
        % 
        %     % from the x values --> structure matrices are generated
        %     % First we need to norm everything to CPN1
        %     % example here for the first param for u1
        %     % x3/(abs(x2)+abs(x3)) * (2*myHeav(x(2)-1)
        %     % using some number examples 
        %     % x2 = 0.8620; x3 = 0.7278
        %     % (0.8620 + 0.7278 * u1) conversion to Norm 1 the following
        %     % (1-abs(0.4578) + 0.4578 * u1) * (1.5899)
        % 
        %     %  (1.5899) here is then used for the phi matrix !
        %     % x1 * (abs(x2) + abs(x3) ......
        %     s(1,1) = obj.darcy_params(3)/(abs(obj.darcy_params(2))+abs(obj.darcy_params(3)))*(2*myHeav(obj.darcy_params(2))-1);
        %     s(2,1) = obj.darcy_params(5)/(abs(obj.darcy_params(4))+abs(obj.darcy_params(5)))*(2*myHeav(obj.darcy_params(4))-1);
        %     s(3,1) = obj.darcy_params(7)/(abs(obj.darcy_params(6))+abs(obj.darcy_params(7)))*(2*myHeav(obj.darcy_params(6))-1);
        %     s(1,2) = obj.darcy_params(10)/(abs(obj.darcy_params(9))+abs(obj.darcy_params(10)))*(2*myHeav(obj.darcy_params(9))-1);
        %     s(2,2) = obj.darcy_params(12)/(abs(obj.darcy_params(11))+abs(obj.darcy_params(12)))*(2*myHeav(obj.darcy_params(11))-1);
        %     %s(:,3) = 0;
        %     phi(1,1) = obj.darcy_params(1)*(abs(obj.darcy_params(2))+abs(obj.darcy_params(3)))*(abs(obj.darcy_params(4))+abs(obj.darcy_params(5)))*(abs(obj.darcy_params(6))+abs(obj.darcy_params(7)))*(2*myHeav(obj.darcy_params(2))-1)*(2*myHeav(obj.darcy_params(4))-1)*(2*myHeav(obj.darcy_params(6))-1);
        %     phi(1,2) = obj.darcy_params(8)*(abs(obj.darcy_params(9))+abs(obj.darcy_params(10)))*(abs(obj.darcy_params(11))+abs(obj.darcy_params(12)))*(2*myHeav(obj.darcy_params(9))-1)*(2*myHeav(obj.darcy_params(11))-1);
        %     %phi(1,3) = 1e-2;
        %     % Now generate the hyCPN1
        %     H = hyCPN1();
        %     % To simplify put everything as an input.
        %     % It is not so problematic as it will be connected to other
        %     % models on the next steps
        %     H.F.input = [s(1:3,:)];
        %     %H.F.algebraic = [s(2:end,:)];
        %     H.phi.equality = [phi];
        %     darcyFactor = dmss(H,0);
        %     darcyFactor.inputName = ["mass_flow_abs","mass_flow_square","turb_friction_factor"];
        % end

        function frictionFactor = frictionFactor(~)
            % Calculate the friction factor for the whole flow region
            % Contains both laminar and turbulent region
            % Using the boolean variable a "toggle" can be made here
            % z1*y1+(1-z1)*y2 = y3 --> eq(1)
            % z1 = boolean to differentiate laminar and turbulent flow
            % y1 = laminar friction factor
            % y2 = turbulent friction factor
            % y1    0	0	1	0
            % y2    1	0	0	1
            % y3    0   1   0   0
            % z1    0   0   1   1
            H1 = hyCPN1();
            H1.F.algebraic = [[sparse([2,3,1,2],[1,2,3,4],1,3,4)]];
            %H1.F.input = [[sparse(1,4,1,1,4)]];
            H1.F.boolean = [[sparse([1,1],[3,4],1,1,4)]];
            H1.phi.equality = [1,-1,1,-1];
            frictionFactor = dmss(H1,0);
            frictionFactor.algebraicName = ["lam_friction_factor","turb_friction_factor","friction_factor"];
            frictionFactor.booleanName = ["Re_bool"];
        end

        function massFlowAbs = massFlowAbs(obj)
            % To model both flow directions, a boolean variable denoting
            % the direction of the mass flow has to be generated. 
            % This MTI model deals with this representation.
            % Everything is put as input as this submodel will be connected
            % with other submodel
            % u4*u1-(1-u4)*u1 = -u2 --> eq (1)
            % eq 1 ensures the value of y1 is the absolute value of u1
            % (2*u4-1)*u1 <= 0--> eq (2)
            % eq 2 give the value of z1 (boolean) depending on the value of
            % u1
            % u3=u2 --> eq (3)
            % eq 3 creates a second variable that is equal to first
            % variable
            % u1    0   1   1   0
            % z1    0   0   1   0
            % u2    1   0   0   0
            % u3    0   0   0   1
            H1 = hyCPN1();
            H1.F.input = [[sparse([2,1,1,4,3],[1,2,3,3,4],1,4,4)]];
            %H1.F.boolean = [[sparse(1,3,1,1,4)]];

            H1.phi.equality = [[1, -1, 2, 0; -1, 0, 0, 1]];
            %H1.phi.inequality = [[0, -1, 2, 0]];
            massFlowAbs = dmss(H1, 0);
            massFlowAbs.inputName = [obj.pipeName+"_m_inlet","mass_flow_abs","mass_flow_square","mass_flow_bool"];
        end
        
        function pressureLoss = pressureLoss(obj)
            % Calculation of overall pressure loss depending on the mass
            % flow, friction factor and inflow pressure
            % y3-u1 + factor *y1*y2*y4*(2*u2-1) = 0
            % u1 = inflow pressure
            % u2 = mass flow bool --> For direction of mass flow
            % y1 = mass flow absolute
            % y2 = friction factor
            % y3 = outlow pressure
            % y4 = mass flow absolute
            %{
            factor = 8*obj.l/(obj.rho*pi()*pi()*power(obj.di,5));

            H1 = hyCPN1();
            H1.F.input = [[sparse([1,2,2],[2,3,4],1,2,4)]];
            H1.F.algebraic = [[sparse([2,1,3,1,3],[1,3,3,4,4],1,3,4)]];
            H1.F.boolean = [[sparse([1],[4],1,1,4)]];
            H1.phi.equality = [[1, -1, factor, -2*factor]];
            pressureLoss = dmss(H1, 0);
            % Put the correct names to ensure the connection is possible
            pressureLoss.algebraicName = ["mass_flow_square","pressure_"+obj.node_target,"mass_flow_abs"];
            pressureLoss.booleanName = ["mass_flow_bool"];
            pressureLoss.inputName = ["pressure_"+obj.node_source,"friction_factor"];
            %}
            factor = 8*obj.l/(obj.rho*pi()*pi()*power(obj.di,5));
            H1 = hyCPN1();
            H1.F.input = [[sparse([1,2,2,3,4],[2,3,4,4,1],1,4,4)]];
            H1.F.algebraic = [[sparse([1,2,1,2],[3,3,4,4],1,2,4)]];
            %H1.F.boolean = [[sparse([1],[4],1,1,4)]];
            H1.phi.equality = [[1, -1, factor, -2*factor]];
            pressureLoss = dmss(H1, 0);
            pressureLoss.algebraicName = ["mass_flow_square","mass_flow_abs"];
            %pressureLoss.booleanName = ["mass_flow_bool"];
            pressureLoss.inputName = [obj.pipeName+"_p_inlet","friction_factor","mass_flow_bool",obj.pipeName+"_p_outlet"];
        end


        % function basePressureLoss = basePressureLoss(obj)
        %     % use the following function for the pressure loss calculation of the model
        %     % y1-u1 = 0 --> eq(1)
        %     % y3-u2 + factor * y1 * y2 * u1 = 0  --> eq(2)
        %     % y1 = U1 = mass flow --> Required to define the square of mass
        %     % flow for pressure loss
        %     % u2 = inflow pressure
        %     % y2 = friction factor (darcy)
        %     % y3 = outflow pressure
        %     % 
        %     %  The structure matrix of the hyCPN1 looks like this
        %     %  u1   1	0	0	1	0
        %     %  u2   0	1	0	0	0
        %     %  y1   0	0	1	1	0
        %     %  y2   0	0	0	1	0
        %     %  y3   0	0	0	0	1
        %     %  The phi matrix then looks like this
        %     %  -1   0   1   0   0   --> eq(1)
        %     %   0   -1  0   factor  1  --> eq(2)
        % 
        %     % Start with a hyCPN1 tensor
        %         H1 = hyCPN1();
        %     % To make everything a bit faster, manually generate the sparse
        %     % matrix required
        %     % Ensure 3 rows and 5 columns --> Bottom 3 
        %         H1.F.algebraic = [sparse([1,2,3,1],[3,4,5,4],1,3,5)];
        %     % Ensure 2 rows and 5 columns --> Top 2 
        %         H1.F.input = [sparse([1,2,1],[1,2,4,],1,2,5)];
        %     % Due to its size for now just generate everything manually
        %         H1.phi.equality = [[-1,0,1,0,0;0,-1,0,8*obj.l/(obj.rho*pi*pi*power(obj.di,5)),1]];
        %     % Now that the tensor is ready, create dmss  
        %         basePressureLoss = dmss(H1,0);
        %     % Put the correct names to ensure the connection is possible
        %         basePressureLoss.algebraicName = ["mass_flow_square", "friction_factor","outflow_pressure"];
        %         basePressureLoss.inputName = ["mass_flow", "inflow_pressure"];
        % end

        function TotalPressureModel(obj)
            % Convenient function to call the total pressure loss model of
            % the pipe
            obj.pressureModel = connect(obj.booleanLamTurb(),...
                obj.darcyLaminar(),...
                obj.darcyFactor(),...
                obj.frictionFactor(), ...
                obj.massFlowAbs(), ...
                obj.pressureLoss());
        end

        function TotalTemperatureModel(obj)

            obj.temperatureModel = connect(obj.energyEquations(),...
                obj.temperatureConnector());
        end

        function temperatureConnector=temperatureConnector(obj)
            H1 = hyCPN1();
            H1.F.input = eye(2);
            H1.phi.equality = [1, -1];
            temperatureConnector = dmss(H1,0);
            temperatureConnector.inputName = [obj.pipeName+"_t_"+obj.N , obj.pipeName+"_t_outlet"];
        end

        function energyEquations=energyEquations(obj)
            % Disaggregated pipe temperature model
            % Minimum disaggregation (N) amount is 2
            % Takes advantage of the "sparsity" of the matrices
            % The matrices are generated as "sparse"
            % Example of structure matrix of a pipe with 2 segments 
            % x --> Temperature (state)
            % xp --> Temperature (stateDerivative)
            % y1 --> Outflow temperature
            % u1 --> inflow mass
            % u2 --> inflow temperature
            % xp1   1	0	0	0	0	0	0	0	0
            % xp2   0	1	0	0	0	0	0	0	0
            % x1    0	0	1	0	0	0	1	0	0
            % x2    0	0	0	1	0	0	0	1	0
            % y1    0	0	0	0	1	0	0	0	0
            % u1    0	0	0	0	0	1	1	1	0
            % u2    0	0	0	0	0	1	0	0	0
            % Always makes use of the ordering of 
            % stateDerivative
            % state
            % algebraic (output)
            % input
            % Example of phi matrix
            % Cth   0	dx*k	0	    0	-cp	cp  0	-Tout*dx*k
            % 0	    Cth	0	    dx*k	0	0	-cp	cp	-Tout*dx*k
            % 0	    0	0	    -1	    1	0	0	0	0   --> (3)
            % Equation (1) and (2) are the energy equation for each segment
            % makes use of the outside temperature as constant
            % If more pipe segments are requiered then equation 1 and 2 are
            % extended further. the "diagonal" parts are continued
            % Equation 3 is always constant
            % (3) connect the last temeperature to an output
            
            % % Ensure the amount of Disaggregation
            % % N < 2 then force it to be 2
            % obj.N = floor(obj.l/obj.delta_x);
            % if obj.N < 2, obj.N=2; end
            % 
            % % Calculate the length of each segment
            % obj.delta_x = obj.l/obj.N;

            % Thermal capacitance of the fluid (water). 
            % depending on the length (segment) and diameter of the pipe
            C_th =  pi()*obj.di*obj.di/4*obj.rho*obj.cp*obj.delta_x; % J/K - Thermal capacitance
            if obj.usePipeWall, C_th = C_th + pi()*((obj.di+2*obj.pipeThickness).^2 - obj.di.^2)/4*obj.rhoPipe*obj.cpPipe*obj.delta_x; end
            
            % Structure matrix for the hyCPN1.
            % First makes use of the "eye" (identity matrix) structure
            % Then add the remaining other parts
            % 
            jA = [1:obj.N*2+1, ... This part is the "diagonal" matrix for xd1 up to u1
                obj.N*2+1, ... This part is for the u2 
                obj.N*2+2:obj.N*3+1,... This part is for the "second" diagonal of x1 to xN
                obj.N*2+2:obj.N*3+1 ... This part is for the horizontal of u1
                ];
            iA = [1:obj.N*2+1, ... This part is the "diagonal" matrix for xd1 up to u1
                obj.N*2+2, ... This part is for the u2 
                obj.N+1:obj.N*2 ...This part is for the "second" diagonal of x1 to xN
                (obj.N*2+1)*ones(1,obj.N) ... This part is for the horizontal of u1
                ];
            % Now generate sparse matrix based on the rows and cols
            % The values are all 1
            % The size of the sparse is "fixed" based on the N segments
            smat = sparse(iA,jA,1,obj.N*2+2,obj.N*3+2);
            
            % Now lets make the component of the phi matrix step by step
            % The locations of ones in the matrix is always consistent
            % iPhi_1 = [obj.N+1];
            % jPhi_1 = [(obj.N*2)+1];
            % vPhi_1 = [1];
            
            % The thermal capacitance parts are also consistent
            % This is a "mini" diagonal depending on the amount of 
            % the segments
            %iPhi_cth = [1:obj.N];
            multiplicator = ones(1,obj.N);

            jPhi_cth = 1:obj.N;
            vPhi_cth = C_th*multiplicator;
            
            % dx*k are also a diagonal
            %iPhi_dxk = [1:obj.N];
            jPhi_dxk = obj.N+1:obj.N*2;
            vPhi_dxk = obj.delta_x*obj.U*multiplicator;
            
            % this is the negative cp component of the equation
            %iPhi_neg_cp = [1:obj.N];
            jPhi_neg_cp = obj.N*2+1:obj.N*3;
            vPhi_neg_cp = -obj.cp*multiplicator;

            % this is the positive cp component of the equation
            %iPhi_cp = [1:obj.N];
            jPhi_cp = obj.N*2+2:obj.N*3+1;
            vPhi_cp = obj.cp*multiplicator;
            
            % This is the constant
            % Always to be placed on the end 
            %iPhi_const = [1:obj.N];
            jPhi_const = (obj.N*3+2)*multiplicator;
            vPhi_const = (-obj.Tout*obj.delta_x*obj.U)*multiplicator;
            
            % Now everything can be combined together
            % This can be generated this way
            % iPhi = [obj.N+1, iPhi_1, iPhi_cth, iPhi_dxk, iPhi_neg_cp, iPhi_cp, iPhi_const];
            % Because of the structure --> iPhis are always 1:obj.N 
            % This can be further simplified
            iPhi = [repmat(1:obj.N,1,5)];
            jPhi = [ jPhi_cth, jPhi_dxk,jPhi_neg_cp, jPhi_cp,jPhi_const];
            vPhi = [vPhi_cth, vPhi_dxk,vPhi_neg_cp, vPhi_cp,vPhi_const];
            
            % Generate the phiMat and put it to the tensor object
            % can be skipped and put directly to the CPN tensor
            %phiMat = sparse(iPhi,jPhi,vPhi, obj.N+2, obj.N*2+4+obj.N+1);
            
            % Generate empty CPN1
            H1 = hyCPN1();
            % Put the "correct" parts of the tensor 
            H1.F.input = smat(obj.N*2+1:obj.N*2+2,:);
            H1.F.stateDerivative = smat(1:obj.N,:);
            H1.F.state = smat(obj.N+1:obj.N*2,:);
            % Put the sparse matrix here directly 
            %H1.phi.equality = phiMat; 
            H1.phi.equality = sparse(iPhi,jPhi,vPhi, obj.N, obj.N*3+2);

            energyEquations = dmss(H1,0);
            energyEquations.inputName = [obj.pipeName+"_m_inlet",obj.pipeName+"_t_inlet"];
            energyEquations.stateName = obj.pipeName+"_t_" + string(1:obj.N);
        end

        function massEquation = massEquation(obj)
            % Generate mass equation
            H1 = hyCPN1();
            H1.F.input = [0,1;1,0];
            H1.phi.equality = [1,1];
            massEquation = dmss(H1,0);
            massEquation.inputName = [obj.pipeName+"_m_inlet", obj.pipeName+"_m_outlet"];
            obj.massModel = massEquation;
        end
    end
end
%     end
% 
%         end
% 
%         function basePressureLoss = BPL(factor)
%             % use the following function for the pressure loss calculation of the model
%             % y1-u1 = 0 --> eq(1)
%             % y3-u2 - factor * y1 * y2 * u1 = 0  --> eq(2)
%             % y1 = U1 = mass flow --> Required to define the square of mass
%             % flow for pressure loss
%             % u2 = inflow pressure
%             % y2 = friction factor (darcy)
%             % y3 = outflow pressure
%             % 
%             %  The structure matrix of the hyCPN1 looks like this
%             %  u1   1	0	0	1	0
%             %  u2   0	1	0	0	0
%             %  y1   0	0	1	1	0
%             %  y2   0	0	0	1	0
%             %  y3   0	0	0	0	1
%             %  The phi matrix then looks like this
%             %  -1   0   1   0   0   --> eq(1)
%             %   0   -1  0   factor  1  --> eq(2)
% 
%             % Start with a hyCPN1 tensor
%                 H1 = hyCPN1();
%             % To make everything a bit faster, manually generate the sparse
%             % matrix required
%             % Ensure 3 rows and 5 columns --> Bottom 3 
%                 H1.F.algebraic = [sparse([1,2,3,1],[3,4,5,4],1,3,5)];
%             % Ensure 2 rows and 5 columns --> Top 2 
%                 H1.F.input = [sparse([1,2,1],[1,2,4,],1,2,5)];
%             % Due to its size for now just generate everything manually
%                 H1.phi.equality = [[-1,0,1,0,0;0,-1,0,sym("factor"),1]];
%             % Now that the tensor is ready, create dmss  
%                 pressureLossModel = dmss(H1,0);
%             % Put the correct names to ensure the connection is possible
%                 pressureLossModel.algebraicName = ["mass_flow_square", "friction_factor","outflow_pressure"];
%                 pressureLossModel.inputName = ["mass_flow", "inflow_pressure"];
%         end
% 
%         function  pressureLossSys = pressureLoss(x, rho, di, l)
%             myHeav = @(X)X>=0;
%             s(1,1) = x(3)/(abs(x(2))+abs(x(3)))*(2*myHeav(x(2))-1);
%             s(2,1) = x(5)/(abs(x(4))+abs(x(5)))*(2*myHeav(x(4))-1);
%             s(3,1) = x(7)/(abs(x(6))+abs(x(7)))*(2*myHeav(x(6))-1);
%             s(1,2) = x(10)/(abs(x(9))+abs(x(10)))*(2*myHeav(x(9))-1);
%             s(2,2) = x(12)/(abs(x(11))+abs(x(12)))*(2*myHeav(x(11))-1);
%             phi(1,1) = x(1)*(abs(x(2))+abs(x(3)))*(abs(x(4))+abs(x(5)))*(abs(x(6))+abs(x(7)))*(2*myHeav(x(2))-1)*(2*myHeav(x(4))-1)*(2*myHeav(x(6))-1);
%             phi(1,2) = x(8)*(abs(x(9))+abs(x(10)))*(abs(x(11))+abs(x(12)))*(2*myHeav(x(9))-1)*(2*myHeav(x(11))-1);
% 
%             H = hyCPN1();
%             H.F.input = [s(1:3,:)];
%             %H.F.algebraic = [s(3,:)];
%             H.phi.equality = [phi];
%             darcyFactor = dmss(H,0);
%             %vpa(darcyFactor.symbolicEquations)
%             darcyFactor.algebraicName = [];
%             darcyFactor.inputName = ["mass_flow","mass_flow_square","friction_factor"];
% 
%             eq = {"y1 = u1";...
%                 "y3 = factor*y1*y2*u1";...
%                 "y4 = u2 - y3"};
%             symbolicParameters = [sym("factor")];
%             numericParameters = [8*l/(rho*pi*pi*power(di,5))];
%             pressureLossModel = sym2dmss(eq, 0);
%             pressureLossModel.algebraicName = ["mass_flow_square", "friction_factor","pressure_loss", "outflow_pressure"];
%             pressureLossModel.inputName = ["mass_flow", "inflow_pressure"];
%             pressureLossModel = pressureLossModel.replaceSymbolicParameters(symbolicParameters, numericParameters);
% 
%             pressureLossSys = connect(darcyFactor,pressureLossModel);
%         end
%     end
% end