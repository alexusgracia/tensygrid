classdef dataProcessing
    % 
    % Utilities for data processing
    %   
    
    % Enrico Uhlenberg, Leona Schnelle - 12/06/2024
    % Adapted by Sarah Casura - February 2025
    methods (Static)
        [stateTrajectoryScaled, inputTrajectoryScaled, ...
          transformation, offset] = ...
            scaleData(stateTrajectory, inputTrajectory, ...
                      stateLowerBound, stateUpperBound, ...
                      inputLowerBound, inputUpperBound)
        [scaleFactor, offset] = ...
            calculateScaleOffset(trajectory, lowerBound, upperBound)
        [H] = heavisideStepFunction(x)
    end
end