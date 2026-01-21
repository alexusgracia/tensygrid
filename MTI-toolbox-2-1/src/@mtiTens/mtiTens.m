classdef (Abstract) mtiTens < handle & matlab.mixin.Heterogeneous
    % Enrico Uhlenberg - 12/06/2024

    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = CPN1();
        end
    end
    
    methods (Abstract)
        result = processXU(sys,xu)
    end


end