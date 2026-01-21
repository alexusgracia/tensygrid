classdef hyCPN1
    %hyCPN1 Stores the model paramters of a dmss model in the structure of a hybrid norm-1 CP
    % decomposed tensor.
    %      H = hyCPN1() Create empty hyCPN1 tensor with hybrid parameter structure.
    %
    % For detailed documentation see <a href="matlab:open((which('hyCPN1Doc.html')))">here</a>

    % Torben Warnecke - 11/06/2024

    properties
        F
        phi
        R
    end
    properties (Hidden)
        sparseIndices
        %         HFxu
        %         HFdxyz
        %         HphiEq
    end

    methods
        function H = hyCPN1()
            %hyCPN1 Stores the model paramters of a dmss model in the structure of a hybrid norm-1 CP
            % decomposed tensor.
            %      H = hyCPN1() Create empty hyCPN1 tensor with hybrid parameter structure.
            %
            % For detailed documentation see <a href="matlab:open((which('hyCPN1Doc.html')))">here</a>

            % Torben Warnecke - 11/06/2024
            
            H.F.state = [];
            H.F.stateDerivative = [];
            H.F.algebraic = [];
            H.F.boolean = [];
            H.F.input = [];

            H.phi.equality = [];
            H.phi.inequality = [];

            H.R = 0;
        end

        function H = set.F(H, newVal)
            H.F = newVal;
            R = [];
            if isfield(H.F, 'state')
                H.F.state = H.trueFalseStructure(H.F.state);
                H.F.state = H.sparseTables(H.F.state);
                H.sparseIndices.state = H.getSparseIndices(H.F.state);
                R = [R; H.sparseIndices.state.nCols];
            end
            if isfield(H.F, 'stateDerivative')
                H.F.stateDerivative = H.trueFalseStructure(H.F.stateDerivative);
                H.F.stateDerivative = H.sparseTables(H.F.stateDerivative);
                H.sparseIndices.stateDerivative = H.getSparseIndices(H.F.stateDerivative);
                R = [R; H.sparseIndices.stateDerivative.nCols];
            end
            if isfield(H.F, 'input')
                H.F.input = H.trueFalseStructure(H.F.input);
                H.F.input = H.sparseTables(H.F.input);
                H.sparseIndices.input = H.getSparseIndices(H.F.input);
                R = [R; H.sparseIndices.input.nCols];
            end
            if isfield(H.F, 'algebraic')
                H.F.algebraic = H.trueFalseStructure(H.F.algebraic);
                H.F.algebraic = H.sparseTables(H.F.algebraic);
                H.sparseIndices.algebraic = H.getSparseIndices(H.F.algebraic);
                R = [R; H.sparseIndices.algebraic.nCols];
            end
            if isfield(H.F, 'boolean')
                H.F.boolean = H.trueFalseStructure(H.F.boolean);
                H.F.boolean = H.sparseTables(H.F.boolean);
                H.sparseIndices.boolean = H.getSparseIndices(H.F.boolean);
                R = [R; H.sparseIndices.boolean.nCols];
            end
            H.R = max(R);%H.sparseIndices.stateDerivative.nCols, H.sparseIndices.state.nCols, H.sparseIndices.input.nCols, H.sparseIndices.algebraic.nCols, H.sparseIndices.boolean.nCols);
        end

        function H = set.phi(H, newVal)
            H.phi = newVal;
            if isfield(H.phi, 'equality')
                H.phi.equality = H.trueFalseParameter(H.phi.equality);
                H.phi.equality = H.sparseTables(H.phi.equality);
                H.sparseIndices.equality = H.getSparseIndices(H.phi.equality);
            end
            if isfield(H.phi, 'inequality')
                H.phi.inequality = H.trueFalseParameter(H.phi.inequality);
                H.phi.inequality = H.sparseTables(H.phi.inequality);
                H.sparseIndices.inequality = H.getSparseIndices(H.phi.inequality);
            end
        end

        function Min = trueFalseStructure(~,Min)
            if isfield(Min, 't') ~= 1
                Mold = Min;
                Min = [];
                if isstring(Mold)
                    Min.t = ((Mold=="true") + (Mold=="1"))==1;
                    Min.f = Mold=="false";
                    Min.c = replace(Mold, ["1", "true", "false"], "0");
                    Min.c(Min.c=="") = "0";
                    if sum(isnan(str2double(Min.c))) == 0
                        Min.c = str2double(Min.c);
                    end
                elseif class(Mold) == "sym"
                    Min.t = zeros(size(Mold));
                    Min.f = zeros(size(Mold));
                    Min.c = Mold;
                    if sum(isSymType(Min.c,'number'),'all')>0
                        Min.c(isSymType(Min.c,'number')) = double(Min.c(isSymType(Min.c,'number')));
                    end
                elseif isnumeric(Mold)
                    Min.t = Mold==1;
                    Min.f = zeros(size(Mold));
                    Min.c = Mold- Min.t;
                else
                    error("Unsupported Structure Matrix type. Must be numeric (e.g. double), string or symbolic.")
                end
            end
        end

        function [Min] = trueFalseParameter(~,Min)
            if isfield(Min, 't') ~= 1
                Mold = Min;
                Min = [];
                if isstring(Mold)
                    Min.t = ((Mold=="true") + (Mold=="1"))==1;
                    Min.f = ((Mold=="false") + (Mold=="-1"))==1;
                    Min.c = replace(Mold, ["1", "-1", "true", "false"], "0");
                    Min.c(Min.c=="") = "0";
                    if sum(isnan(str2double(Min.c))) == 0
                        Min.c = str2double(Min.c);
                    end
                elseif class(Mold) == "sym"
                    Min.t = zeros(size(Mold));
                    Min.f = zeros(size(Mold));
                    Min.c = Mold;
                    if sum(isSymType(Min.c,'number'),'all')>0
                        Min.c(isSymType(Min.c,'number')) = double(Min.c(isSymType(Min.c,'number')));
                    end
                elseif isnumeric(Mold)
                    Min.t = Mold==1;
                    Min.f = Mold==-1;
                    Min.c = Mold- Min.t - Min.f*(-1);
                else
                    error("Unsupported Parameter Matrix type. Must be numeric (e.g. double), string or symbolic.")
                end
            end
        end

        function Min = sparseTables(~,Min)
            if issparse(Min.t) ~= 1
                %Mold = Min;
                Min.t = sparse(Min.t);
            end
            if issparse(Min.f) ~= 1
                Min.f = sparse(Min.f);
            end
            if class(Min.c) ~= "sym" && issparse(Min.c)~=1
                Min.c = sparse(Min.c);
            else
                Min.c = Min.c;
            end
        end
    end
end

