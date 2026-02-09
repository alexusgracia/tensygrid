classdef otvlTens < mtiTens
    % OTVLTENS Construct a otvlTens object.
    % Multidimensional structure containing the state equations
    % of a multilinear model which is representable by OTVLs. 
    %
    % sys = otvlTens(FPhi, Fa, Fc, varargin) constructs an otvlTens
    % object where varargin contains objects of type otvl, which 
    % represent the equation structure of the multilinear state 
    % equations.
    %
    % inputs: 
    %   FPhi: scaling parameter for otvl addends
    %   Fa: shifting parameter for individual variables
    %   Fc: offset parameter for individual otvls
    %   varargin: variable number of otvl objects
    %
    % properties:
    %   structure: cell array containing otvls
    %   FPhi: scaling parameter for otvl addends
    %   Fa: shifting parameter for individual variables
    %   Fc: offset parameter for individual otvls
    %   equationCount: number of state equations
    %   inputCount: number of inputs (currently fixed to 0)
    %
    % For detailed documentation see <a
    % href="matlab:open((which('otvlTensDoc.html')))">here</a>.
    %
    % See also otvl mlgreyest

    % Marah Engels - 11/06/2024

    properties
        structure %contains TVLs
        FPhi %scaling parameter for TVL addends
        Fa %shifting parameter referring to individual variables
        Fc %Offset Parameter for individual TVLs
       
    end

    properties (Access = protected)
        NVariables %number of variables in each TVL (might be obsolete)
        NEquations %number of equations
        NInputs %number of inputs
        NRows %number of rows per TVL
        NRowsTotal %total number of rows
        DontCareArray %2D array of TVL Don'tCare-dimension
        BooleanArray %2D array of TVL Boolean-dimension
        IndexVector %index vector allocating rows of DontCareArray and Boolean Array to TVLs
        FaInter %internal shifting parameter 
    end

    properties (Dependent, SetAccess = private)
        equationCount int64
        inputCount int64
    end

    methods 
        %constructor method
        function sys = otvlTens(FPhi, Fa, Fc, varargin)
            % OTVLTENS Construct a otvlTens object.
            % Multidimensional structure containing the state equations
            % of a multilinear model which is representable by OTVLs. 
            %
            % sys = otvlTens(FPhi, Fa, Fc, varargin) constructs an otvlTens
            % object where varargin contains objects of type otvl, which 
            % represent the equation structure of the multilinear state 
            % equations.
            %
            % inputs: 
            %   FPhi: scaling parameter for otvl addends
            %   Fa: shifting parameter for individual variables
            %   Fc: offset parameter for individual otvls
            %   varargin: variable number of otvl objects
            %
            % properties:
            %   structure: cell array containing otvls
            %   FPhi: scaling parameter for otvl addends
            %   Fa: shifting parameter for individual variables
            %   Fc: offset parameter for individual otvls
            %   equationCount: number of state equations
            %   inputCount: number of inputs (currently fixed to 0)
            %
            % For detailed documentation see <a
            % href="matlab:open((which('otvlTensDoc.html')))">here</a>.
            %
            % Author(s): Marah Engels
            % Copyright and Disclaimers: see mti.systems (2024)
            %
            % See also otvl mlgreyest
            
            if nargin < 1 %create object with empty properties
                    sys.structure = {};
                    sys.FPhi = [];
                    sys.Fa = [];
                    sys.FaInter = [];
                    sys.Fc = [];
                    sys.NVariables = [];
                    sys.NEquations = [];
                    sys.NRows = [];
                    sys.NRowsTotal = [];
                    sys.DontCareArray = [];
                    sys.BooleanArray = [];
                    sys.IndexVector = [];
                    
            elseif nargin < 5
                error('Number of inputs not sufficient.')
            else 
                %check type of given parameters
                if ~(isa(FPhi, 'double')  && isa(Fa, 'double') && isa(Fc, 'double'))
                    error(['Parameters FPhi, Fc and Fa have to be given as an array. ' ...
                        'If you want to set all Parameters to a value of 1, ' ...
                        'call otvlmss([],[],[],[], varargin'])
                else 
                    %initialize variables
                    otvlElements = varargin;
                    nElements= size(otvlElements,2);   
                    FPhi_temp = cell(nElements,1);
                    sys.NVariables = zeros(nElements,1);
                    sys.NRows = zeros(nElements,1);

                    %read information from given OTVLs
                    for iElement = 1:1:nElements
                        element = otvlElements{iElement};
                        if ~isa(element,'otvl') %check data type
                            typeException = MException('otvlmss:WrongInputType','Wrong input type at position %i. All variable input elements have to be of class otvl', iElement+4);
                            throw(typeException)
                        else
                            nVariables = size(element.ttable,2); %get number of variables in respective TVL
                            nRows = size(element.ttable,1); %get number of rows in respective TVL
                            FPhi_temp{iElement,1} = ones(nRows,1);
                            sys.NVariables(iElement,1) = nVariables;
                            sys.NRows(iElement,1) = nRows;
                        end                        
                    end

                    if ~(nElements == max(sys.NVariables))
                        error('Number of given OTVL has to equal number of variables')
                    end
                    
                    %assign values to properties
                    sys.structure = otvlElements;
                    sys.FPhi = cell2mat(FPhi_temp);
                    sys.Fa = zeros(2,max(sys.NVariables));
                    sys.FaInter = sys.Fa;
                    sys.FaInter(2,:) = 1-sys.Fa(2,:);
                    sys.Fc = zeros(max(sys.NVariables),1);

                    % inherit parameters from input, if given
                    %TODO: Check whether input parameters have correct
                    %dimension
                    if ~(isempty(FPhi))
                        sys.FPhi = FPhi;
                    end

                    if ~(isempty(Fa))
                        sys.Fa = Fa;
                        sys.FaInter = sys.Fa;
                        sys.FaInter(2,:) = 1-sys.Fa(2,:);
                    end  

                    if ~(isempty(Fc))
                        sys.Fc = Fc;
                    end  
                    
                    %initialize 2D-arrays for simulation and identification of
                    %otvlmss-objects
                    sys.NEquations = size(sys.structure,2);
                    sys.NInputs = 0;
                    sys.NRowsTotal = sum(sys.NRows); % sum of rows in all TVLs of the system
                    sys.IndexVector = zeros(sys.NRowsTotal,1); %index array to allocate TVL rows to equations
                    sys.DontCareArray = false(sys.NRowsTotal, sys.NEquations); % merged Dont-Care-arrays of all TVLs
                    sys.BooleanArray = false(sys.NRowsTotal, sys.NEquations); % merged Boolean-arrays of all TVLs
                    
                    %assign values to 2D-arrays
                    iRow = 1;
                    for iEquation = 1:1:sys.NEquations
                        sys.DontCareArray(iRow:iRow+sys.NRows(iEquation)-1,:) = sys.structure{iEquation}.ttable(:,:,1);
                        sys.BooleanArray(iRow:iRow+sys.NRows(iEquation)-1,:) = sys.structure{iEquation}.ttable(:,:,2);
                        sys.IndexVector(iRow:iRow+sys.NRows(iEquation)-1,:) = iEquation*ones(sys.NRows(iEquation),1);
                        iRow = iRow + sys.NRows(iEquation); 
                    end

                end
            end
            
        end
        
        function result = processXU(obj,xu)
            xPre = (obj.FPhi.*prod(obj.DontCareArray + (~obj.DontCareArray).*(obj.BooleanArray.*(xu' + obj.FaInter(1,:)) ...
                    +(~obj.BooleanArray).*(obj.FaInter(2,:) - xu')),2))'; %% TODO check xu orientation
            result = accumarray(obj.IndexVector, xPre')' + obj.Fc';
        end

        function updateSparseIndicesAndSize(obj)
            obj.structure = obj.structure;
            %placeholder until this is refactored
        end

        function equationCount = get.equationCount(obj)
            equationCount = obj.NEquations;
        end

        function inputCount = get.inputCount(obj)
            inputCount = 0*obj.NInputs;
            %warning('Warning: inputCount is fixed to zero, as only autonomous systems are implented currently');
        end
        %method simulation of otvlmss object
        x = otvlsim(otvl, u, t, x0) %currently two different versions for equation with and without additional parameters
        x = otvlsimA(otvl, u, t, x0)

        %method ALS-algorithm for parameter identification of otvlmss
        %object
        %[Fphi, Falpha, cost] = otvl2Als(data, otvlmssObj, options) %currently two different versions for equation with and without additional parameters
        [Fphi, Flambda, Fa, Fc, cost] = otvl2AlsA(data, otvlmssObj, options)
    end
end