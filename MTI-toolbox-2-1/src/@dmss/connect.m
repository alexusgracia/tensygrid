function sys = connect(sys, varargin)
% CONNECT - connect signals of dmss objects
%
% possible input scenarios:
% sys = connect(sys): 
%   connects an appended system based on signal names
%
% sys = connect(sys, stateConnections, algebraicConnections, booleanConnections, similarInputs):
%   connects an appended system based on indexes (state-input-/
%   algebraic-inputs-/ boolean-input-connections and similar
%   inputs). Indexes are parsed as connection matrices with input index as
%   first value in a row and the variable indices as following values [inputIdx, variableIdx, ...; ...].
%   If more then 2 values are in a row a additional (linear) connection equation will be added: u = var(1) + var(2) + ...
%   If a variable index is parsed as negative value it is subtracted in the
%   connection equation: u = -var(1) ...
%
% sys = connect(sys, sys2, ..., sysN):
%   connects multiple systems based on signal names
%
% For detailed documentation see <a href="matlab:open((which('connectDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

algebraicConnections = [];
booleanConnections = [];
similarInputs = [];

if nargin == 1
% Connect appended system by names
    sys = sys.connectByNames;

elseif nargin > 1
    if isnumeric(varargin{1,1})
% Connect appended system by indices
        stateConnections = varargin{1,1};
        if nargin > 2
            algebraicConnections = varargin{1,2};
        end
        if nargin > 3
            booleanConnections = varargin{1,3};
        end
        if nargin == 5
            similarInputs = varargin{1,4};
        end
        if nargin > 5
            error('Too many input arguments.')
        end
        
        sys = sys.connectByIndices(stateConnections, algebraicConnections, booleanConnections, similarInputs);

    elseif class(varargin{1,1}) == "dmss"
        % Connect multiple systems by names
        sys = append(sys, varargin{:});
        sys = sys.connectByNames;
    else
        error('Wrong input data type. Inputs must be dmss-class object or integer (connection) matrices.')
    end
end

sys = sys.trivialReduction;

end

