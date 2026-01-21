
classdef mss < mti
    properties (Access = private)
        data  (1,:) {mustBeA(data,"mtiTens")} = [CPN1()] % Array unfortunately needs non-abstract default value.
    end
    properties (Dependent, SetAccess = private)
        F mtiTens
        G mtiTens
        n int64
        m int64
        p int64
        ntype char

    end
    properties (Access = public)
        ts double {mustBeNonnegative}
        stateName string
        inputName string
        outputName string
    end
    
    methods
        function sys = mss(varargin)
            % MSS Construct a mss object
            %   sys = mss(sys) Copy constructor, when sys is a mss object
            %
            %   sys = mss(Ftens) Constructs a mss object from a tensor
            %   given as a CPN1 object. 
            %
            %   sys = mss(Ftens,Gtens) constructs a
            %   mss object from two CPN1 object, with the first
            %   tensor representing state-transition and the second output
            %   transition
            %   
            %   sys = mss(Fmat) Constructs a mss object from a matrix
            %   representing the . 
            %
            %   sys = mss(Fmat,Gmat) constructs a
            %   mss object from two matrices representing a MTI system, with the first
            %   matrix representing state-transition parameters and the
            %   second output-transition-parameters
            %
            %   sys = mss(ss) Constructs a mss based on a given linear
            %   model in the matlab-native ss format. 
            %
            %   sys = mss(otvlTens) Constucts a mss based on a given
            %   ternary vector list. For details see TODO
            %
            %   sys = mss(..., ts) constructs the mss object with a given
            %   timestep size ts, which can be 0 (continuous-time) or bigger
            %   than 0 (discrete-time). ts is always the last (scalar) input
            %   argument of mss. default: ts = 0
            %    
            %   Note: The CPN1 objects can also be passed as ktensor from
            %   "Tensor Toolbox".
            %
            %   For detailed documentation see <a href="matlab:open((which('mssDoc.html')))">here</a>
            %   See also dmss msim
            
            % Enrico Uhlenberg , Torben Warnecke, Carlos Cateriano
            % Yáñez, Christoph Kaufmann, Marah Engels, Leandro Samaniego,
            % Leona Schnelle - 11/06/2024          

            % Environment check
            environmentChecker(false)

            % Empty Constructor
            if(nargin == 0)
                return
            end
            
            % Copy Contructor
            if (nargin == 1 && isa(varargin{1},'mss'))
                sys = copy(varargin{1});
                return
            end

            if nargin > 3
                error('Too many input arguments. Inputs can be F, G and ts!')
            end
            
            n_input_tensors = nargin;

            % set timestep size
            if isnumeric(varargin{end}) && isscalar(varargin{end})
                sys.ts = varargin{end};
                n_input_tensors = n_input_tensors -1; 
            else
                sys.ts = 0;
            end  
            
            if(sys.ts == 0 && isa(varargin{1},'otvlTens'))
                error('Continous otvl systems currently not supported')
            end

            for i = 1:n_input_tensors
                switch(checkInputType(varargin{i}))
                    case 'matrix'
                        try
                            sys.data(i) = CPN1(varargin{i});
                        catch E
                            rethrow(E);
                        end 
                    case 'fullTens'
                        try
                            sys.data(i) = CPN1(varargin{i});
                        catch E
                            rethrow(E);
                        end
                    case 'CPN1'
                        try
                            sys.data(i) = varargin{i};
                        catch E
                            rethrow(E);
                        end
                    case 'ss'
                        [tempF, tempG] = mss.ss2Mss(varargin{i},'1');
                        try
                            sys.data(i) = CPN1(tempF.U,tempF.phi);
                            sys.data(i+1) = CPN1(tempG.U,tempG.phi); % Bit ugly but since the 'ss' is asingle input representing two tensors we have to "skip" forward in the loop here
                        catch E
                            rethrow(E);
                        end
                    case 'ktensor'
                        try
                            sys.data(i) = CPN1(varargin{i});
                        catch E
                            rethrow(E);
                        end
                    case 'otvlTens'
                        try
                            if i > 1
                                error('Currently only Systems with without output equation are supported'); % Enrico: Temporary solution until the concept of a complete Xvl System is understood
                            end
                            sys.data(i) = varargin{i};
                        catch E
                            rethrow(E);
                        end
                    otherwise
                        error('Input type not supported')
                end
            end

            sys.updateSparseIndicesAndSize();
        end
        
        % Get / Set Functions
        function val = get.F(sys)
            if(isempty(sys.data))
                val = mtiTens.empty();
                return
            end
            val = sys.data(1);            
        end

        function val = get.G(sys)
            if(length(sys.data) == 1)
                val = mtiTens.empty();
                return
            end
            val = sys.data(2);

        end
            
        function set.F(sys,input)
            sys.data(1)= input;
            sys.updateSparseIndicesAndSize();
        end

        function set.G(sys,input)
            sys.data(2) = input;
            sys.updateSparseIndicesAndSize();
        end
        
        function n = get.n(sys)
            n = sys.data(1).equationCount;
        end

        function m = get.m(sys)
            m = sys.data(1).inputCount;
        end

        function p = get.p(sys)
            if(length(sys.data) == 1)
                    %p = sys.n; % is that right?
                    p = 0;
            else
                    p = sys.data(2).equationCount;
            end
        end
		
        function val = get.ntype(sys)
            % Type Verification of CPN1 
            if (isa(sys.data(1), 'CPN1'))
                val = '1';
                return
            end

%             if (isa(sys.data(1), 'norm2CPN'))
%                 val = '2';
%                 return
%             end
            %Type Verification Extend: p,z -> Not implemented yet
            %if (isa(sys.data(1), 'nameOfz'))
            %   val = 'z';
            %end

            %if (isa(sys.data(1), 'nameOfp'))
            %   val = 'p';
            %end
        end

        function set.stateName(sys,names)
            arguments
                sys
                names (1,:) string
            end
            if length(names) ~= sys.n
                error('State name vector must be a string array with as many names as the model has number of states (sys.n).')
            else
                sys.stateName = names;
            end
        end

        function set.inputName(sys,names)
            arguments
                sys
                names (1,:) string
            end
            if length(names) ~= sys.m
                error('Input name vector must be a string array with as many names as the model has number of inputs(sys.m).')
            else
                sys.inputName = names;
            end
        end

        function set.outputName(sys,names)
            arguments
                sys
                names (1,:) string
            end
            if length(names) ~= sys.p
                error('Output name vector must be a string array with as many names as the model has number of outputs (sys.p).')
            else
                sys.outputName = names;
            end
        end
        % User Exposed Functions

        function newSys = MEXaccelerate(sys,enable)
            arguments
                sys mss
                enable logical 
            end
            if(enable)
                XUfunName = "MTISIM";
            else
                XUfunName = "accumarray";
            end

            sys.F.changeXUfun(XUfunName);
            if(isa(sys.G,"CPN1"))
                sys.G.changeXUfun(XUfunName);
            end
            newSys = sys;
        end
        
        x = getNextState(sys,xu)

        y = getOutput(sys,xu)

        [mtiSystemTransformed] = mss2mss(mtiSystem, transformation, offset)

    end

    methods (Static)
        [F, G] = ss2Mss(lsys,normtype)
        [msys] = mlinearize(model,low_bnd,up_bnd,level,max_order)
        msys = rmss(n,m,p,r)
        msys = drmss(n,m,p,r,Ts)
    end


    % Methods for Simulation
    methods (Access = {?simMTI}) 

        function [tOut, x] = getStateTrajectory(sys, u, t, x0, method)
            if length(t)>1
                tspan = [t(1) t(end)];
            elseif length(t) == 1
                tspan = [0 t(end)];
            end
            
            [tOut, x] = ode45(@(tOut,x)stateHandle(sys,tOut,x,u,t,method), tspan, x0);
            
            function dxdt = stateHandle(sys, tOut, x, u, t, method)
                n = sys.n;
                p = sys.p;
                m = sys.m;
                xu = [];
                for j = 1:n
                    xu = [xu; x(j)];
                end
                for j = 1:m
                    u(j) = interp1(t, u(:,j), tOut, method);
                    xu = [xu; u(j)];
                end
                dxdt = sys.data(1).processXU(xu);
            end

        end

        function updateSparseIndicesAndSize(sys)
            for i=1:numel(sys.data)
                sys.data(i).updateSparseIndicesAndSize();
            end
        end

    end



end

function inputType = checkInputType(unknownInput)
    if (isa(unknownInput,'double') && ~iscell(unknownInput))
        if ndims(unknownInput) <=2
            inputType = 'matrix';
        else
            inputType = 'fullTens';
        end
    elseif isa(unknownInput, 'fullTens')
        inputType = 'fullTens';
    elseif((isa(unknownInput,'CPN1')))
        inputType = 'CPN1';
    elseif(isa(unknownInput,'ktensor'))
        inputType = 'ktensor';
    elseif(isa(unknownInput,'otvlTens'))
        inputType = 'otvlTens';
    elseif(isa(unknownInput,'ss'))
        inputType = 'ss';
    else
        inputType = 'unsupported';
    end   
    
end