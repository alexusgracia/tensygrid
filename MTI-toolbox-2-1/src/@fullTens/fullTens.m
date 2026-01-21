classdef fullTens < mtiTens
    % Enrico Uhlenberg - 12/06/2024
    properties
        Tens double
    end

    methods
        function obj = fullTens(A)
            obj.Tens = A;
        end
    end
end