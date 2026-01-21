classdef dmss
    %DMSS Construct sparse norm-1 implicit/descriptor multilinear
    % time invariant state space model in canonical polyadic decomposition
    % as a dmss object
    %
    %   sys = dmss() Create empty dmss object with an empty hyCPN1
    %   tensor
    %
    %   sys = dmss(H) Construct a continuous-time (ts = 0) dmss
    %   object when H is a hyCPN1 tensor
    %
    %   sys = dmss(mssSys) Convert a mss object into a dmss object
    %   when mssSys is a mss object with a similar stepsize
    %   sys.ts = mssSys.ts
    %
    %   sys = dmss(H,ts) Construct a dmss object when H is a hyCPN1
    %   tensor, with ts as stepsize (ts = 0: continuous-time,
    %   ts > 0: discrete-time)
    %
    % For detailed documentation see <a href="matlab:open((which('dmssDoc.html')))">here</a>

    % Torben Warnecke - 11/06/2024

    properties
        % model tensor ("advanced" sparse)
        H hyCPN1

        % number of variables
        n   % states
        m   % inputs
        p   % algebraic (auxiliaries/outputs)
        q   % boolean
        
        % output indices / Analysis Points
        outputIndex % marks which variables are outputs (e.g. for simulation or connecting)
        
        % analysisPoint

        % number of equalities
        nEq
        
        % number of inequalities
        nIneq
        
        stateName string
        stateUnit string
        
        algebraicName string
        algebraicUnit string
        
        inputName string
        inputUnit string
        
        booleanName string
        booleanUnit string

        % indices of different input types
        manipulatedVariable uint16 
        measuredDisturbance uint16
        unmeasuredDisturbance uint16

        ts double {mustBeNonnegative} % 0 for continuous-time, >0 for discrete-time

        discreteDifferentiation string ="none"
    end

    properties (Hidden)

    end

    methods
        function sys = dmss(H, ts)
            % DMSS Construct an implicit/descriptor MTI model 
            % as a dmss object
            %
            %   sys = dmss() Create empty dmss object with an empty hyCPN1 
            %   tensor
            %
            %   sys = dmss(H) Construct a continuous-time (ts = 0) dmss
            %   object when H is a hyCPN1 tensor
            %
            %   sys = dmss(mssSys) Convert a mss object into a dmss object
            %   when mssSys is a mss object with a similar stepsize
            %   sys.ts = mssSys.ts
            %
            %   sys = dmss(H,ts) Construct a dmss object when H is a hyCPN1
            %   tensor, with ts as stepsize (ts = 0: continuous-time,
            %   ts > 0: discrete-time)
            %
            % For detailed documentation see <a href="matlab:open((which('dmssDoc.html')))">here</a>

            % Torben Warnecke - 11/06/2024



            % Environment check
            environmentChecker(false)
            
            if nargin < 1
                % creates empty object
                %sys.ts = 0;
                sys.H = hyCPN1;

                sys.outputIndex.state = [];
                sys.outputIndex.algebraic = [];
                sys.outputIndex.boolean = [];

            elseif nargin <= 3
                if isa(H, 'mss')
                    % transforms mss-object to dmss-object
                    [sys.H, sys.ts] = sys.mss2dmss(H);

                elseif isa(H, 'hyCPN1')
                    sys.H = H;

                    if nargin >= 2
                        sys.ts = ts;
                    else
                        sys.ts = 0;
                    end
                else
                    error('Data type of H not supported!')
                end

                if nargin >=3
                    if ~isfield(outputIndex, 'state')||~isfield(outputIndex, 'algbraic')||~isfield(outputIndex, 'boolean')
                        error('outputIndex-input must be a struct with the fields state, algebraic and boolean.')
                    else
                        if max(outputIndex.state) > sys.n
                            error('Output inidces for state variables cannot be higher than the number of state variables!')
                        elseif max(outputIndex.algebraic) > sys.p
                            error('Output inidces for algebraic variables cannot be higher than the number of algebraic variables!')
                        elseif max(outputIndex.boolean) > sys.q
                            error('Output inidces for boolean variables cannot be higher than the number of boolean variables!')
                        else
                            %sys.outputIndex = outputIndex;
                        end
                    end
                else
                    sys.outputIndex.state = [];
                    sys.outputIndex.algebraic = [];
                    sys.outputIndex.boolean = [];
                end
                sys.n = max([size(sys.H.F.state.t,1), size(sys.H.F.stateDerivative.t,1)]);
                sys.m = size(sys.H.F.input.t,1);
                sys.p = size(sys.H.F.algebraic.t,1);
                sys.q = size(sys.H.F.boolean.t,1);
                sys.nEq = size(sys.H.phi.equality.t,1);
                sys.nIneq = size(sys.H.phi.inequality.t,1);

                sys = sys.trivialReduction();
            end
        end

        function sys = set.stateName(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.n
                error('State name vector must be a string array with as many names as the model has number of states (sys.n).')
            else
                sys.stateName = strArray;
            end
        end
        function sys = set.stateUnit(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.n
                error('State unit vector must be a string array with as many units as the model has number of states (sys.n).')
            else
                sys.stateUnit = strArray;
            end
        end

        function sys = set.inputName(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.m
                error('Input name vector must be a string array with as many names as the model has number of inputs (sys.m).')
            else
                sys.inputName = strArray;
            end
        end
        function sys = set.inputUnit(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.m
                error('Input unit vector must be a string array with as many units as the model has number of inputs (sys.m).')
            else
                sys.inputUnit = strArray;
            end
        end

        function sys = set.algebraicName(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.p
                error('Algebraic name vector must be a string array with as many names as the model has number of algebraic variables (sys.p).')
            else
                sys.algebraicName = strArray;
            end
        end
        function sys = set.algebraicUnit(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.p
                error('Algebraic unit vector must be a string array with as many units as the model has number of algebraic variables (sys.p).')
            else
                sys.algebraicUnit = strArray;
            end
        end

        function sys = set.booleanName(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.q
                error('Boolean name vector must be a string array with as many names as the model has number of boolean variables (sys.q).')
            else
                sys.booleanName = strArray;
            end
        end
        function sys = set.booleanUnit(sys, strArray)
            arguments
                sys
                strArray (1,:) string
            end
            if length(strArray) ~= sys.q
                error('Boolean unit vector must be a string array with as many units as the model has number of boolean variables (sys.q).')
            else
                sys.booleanUnit = strArray;
            end
        end

    end
    methods (Static, Access = public) 
            [H, ts] = mss2dmss(sys)
            sys = sym2dmss(sys, symfun, Ts)
            sys = poly2dmss(poly_str_cell, Ts)
            sys = str2dmss(str_cell, Ts)
    end
end