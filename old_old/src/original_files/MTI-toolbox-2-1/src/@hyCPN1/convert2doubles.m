function H = convert2doubles(H)
% Torben Warnecke - 11/06/2024
    function [Mout] = doubleMatrix(Min)
        Mout = Min;
        if isfield(Min, 't') == 1
            if (isa(Min.c,'double') == 0)
                try
                    Mout.c = sparse(double(Min.c));
                catch
                    error('Model tensor could not be converted into a doubles.')
                end
            end
        else
            if (isa(Min,'double') == 0)
                try
                    Mout = sparse(double(Min));
                catch
                    error('Model tensor could not be converted into a doubles.')
                end
            end
        end

    end

H.F.stateDerivative = doubleMatrix(H.F.stateDerivative);
H.F.state = doubleMatrix(H.F.state);
H.F.algebraic = doubleMatrix(H.F.algebraic);
H.F.input = doubleMatrix(H.F.input);
H.F.boolean = doubleMatrix(H.F.boolean);

H.phi.equality = doubleMatrix(H.phi.equality);
H.phi.inequality = doubleMatrix(H.phi.inequality);

% if marker == 1
%     disp('Symbolic matrices have been converted to doubles.')
% end
end

