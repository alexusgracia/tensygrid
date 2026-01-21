classdef CPN1 < cpnTens

    properties (Dependent, SetAccess = private)
        equationCount int64
        inputCount int64
    end
    methods
        function obj = CPN1(varargin)
             % CPN1 constructs a norm-1 CP-decomposed tensor object
             %
             %    obj=CPN1(U,phi) takes the structural matrix U and the
             %    parameter matrix phi specified as matrices and returns
             %    CPN1 object. It is assumed that each state has a
             %    state equation. 
             %
             %    obj=CPN1(tens) takes a ktensor from the "Tensor
             %    Toolbox" and returns a CPN1 decomposed tensor.
             %    
             %    obj=CPN1(M) takes a multilinar model represented as 
             %    matrix, vector double array, or double tensor and returns 
             %    a CPN1 object.
             %   
             %    obj=CPN1(A) takes the special case of 2x2x1 tensor
             %    specified by the matrix A and returns a CPN1 object. 
             %
             % Example:   
             %    The first-order multilinear model with one state and one input
             %    xp1=2*x1*u1+3*x1 is created by: 
             %          U=[1 1; 1 0];
             %        phi=[2 3];     
             %        obj=CPN1([1 1; 1 0],[2 3])
             %
             % For detailed documentation see <a href="matlab:open((which('CPN1Doc.html')))">here</a>
             % See also jacobian

             % Enrico Uhlenberg, Carlos Cateriano Yáñez, Christoph Kaufmann, Torben Warnecke, Leona Schnelle - 11/06/2024


            if(nargin == 0)
                obj.U = double.empty();
                obj.phi = double.empty();
            elseif(nargin == 1)
                if (isa(varargin{1}, 'double') && isequal(size(varargin{1}), [2,2]))
                    warning('CPN1:matrix2x2Warning', 'Input array is treated as 2x2-Matrix, if it should be treated as tensor, use the fullTens-class: CPN1(fullTens(A))')
                end

                if isa(varargin{1},'ktensor')
                    [U,phi] = obj.ktensor2CPN1(varargin{1});
                elseif (isa(varargin{1}, 'double') && ndims(varargin{1}) <= 2) % checks if input is a matrix or vector double array
                    [U,phi] = cpnTens.mat2Cpn(varargin{1});
                elseif (isa(varargin{1}, 'double') && ndims(varargin{1}) >2) % checks if input is a double tensor (matlab native)
                    [U,phi] = cpnTens.fulltens2Cpn(varargin{1});
                elseif isa(varargin{1}, 'fullTens') % checks if input is a fullTens-object (for the special case of 2x2x1 Tensors)
                    [U,phi] = cpnTens.fulltens2Cpn(varargin{1}.Tens);
                end


                obj.U = U;
                obj.phi = phi;

                setXUfun(obj);
            elseif(nargin == 2)
                                
                obj.U = varargin{1};
                obj.phi = varargin{2};

                setXUfun(obj);
                % TODO size check
                
            else
                error("Error initializing norm1 Object, invalid number of inputs. Expected parameters U and phi ")
            end

        end

        function result = processXU(sys,xu)
            result = sys.processXUfun(sys,xu);
        end

        function equationCount = get.equationCount(obj)
            equationCount = size(obj.phi,1);
        end

        function inputCount = get.inputCount(obj)
            inputCount = size(obj.U,1) - obj.equationCount;
        end

        function changeXUfun(obj, changeTo)
            arguments
                obj CPN1
                changeTo string
            end
            switch changeTo
                case "accumarray"
                    obj.processXUfun = @(sys,xu) sys.phi*accumarray(sys.ColInd,1+sys.DataVec.*(xu(sys.RowInd)+sys.SignVec),[sys.rowsU 1],@prod,1);
                case "MTISIM"
                    obj.U = sparse(obj.U);
                    obj.phi = sparse(obj.phi);
                    obj.processXUfun = @(sys,xu) MTISIM(sys.U,sys.phi,xu(sys.equationCount+1:end),xu(1:sys.equationCount)); 
                otherwise
                    error("Option not supported by CPN1. Supprted Options: accumarray or MTISIM" )
            end
        end
    end
    methods (Static, Access = public) 
            [U,phi] = ktensor2CPN1(T)
    end

end
function setXUfun(obj)
    % This defines a lambda function as a property of the CPN1 object,
    % circumventing two different subclasses for bool and double
    % calculations. This will probably be obsolete as soon as we progress
    % on sparse/nonsparse differentiation and also when we want to be able
    % to handle hybrid systems
    if(islogical(obj.U))
        obj.processXUfun = @(sys,xu) sys.phi*accumarray(sys.ColInd,xu(sys.RowInd),[sys.rowsU 1],@prod,1);
    else
         obj.processXUfun = @(sys,xu) sys.phi*accumarray(sys.ColInd,1+sys.DataVec.*(xu(sys.RowInd)+sys.SignVec),[sys.rowsU 1],@prod,1);
    end
end


