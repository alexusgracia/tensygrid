function [stateTrajectoryScaled, inputTrajectoryScaled, ...
          transformation, offset] = ...
    scaleData(stateTrajectory, inputTrajectory, ...
              stateLowerBound, stateUpperBound, ...
              inputLowerBound, inputUpperBound)
%
% <scaleData> - scale the state and input data in a range between lower and 
% upper bound; return corresponding transformation matrix and offset vector
%
%
% Input parameters: 
%
%  - stateTrajectory: state trajectory x
%                     This can be a matrix containing several trajectories 
%                     for several states (one trajectory per column).
%
%  - inputTrajectory: input trajectory u
%                     This can be a matrix containing several trajectories 
%                     for several inputs (one trajectory per column).
%
%  - stateLowerBound: lower bound(s) to scale the stateTrajectory to
%                     This can be a scalar (which is applied to all
%                     trajectories in stateTrajectory), or a vector with 
%                     the number of elements corresponding to the number of 
%                     columns in stateTrajectory. 
%
%  - stateUpperBound: upper bound(s) to scale the stateTrajectory to
%                     This can be a scalar (which is applied to all
%                     trajectories in stateTrajectory), or a vector with 
%                     the number of elements corresponding to the number of 
%                     columns in stateTrajectory. 
%
%  - inputLowerBound: lower bound(s) to scale the inputTrajectory to
%                     This can be a scalar (which is applied to all
%                     trajectories in inputTrajectory), or a vector with 
%                     the number of elements corresponding to the number of 
%                     columns in inputTrajectory. 
%
%  - inputUpperBound: upper bound(s) to scale the inputTrajectory to
%                     This can be a scalar (which is applied to all
%                     trajectories in inputTrajectory), or a vector with 
%                     the number of elements corresponding to the number of 
%                     columns in inputTrajectory. 
%
%
% Output parameters: 
%
%  - stateTrajectoryScaled: the state trajectory scaled between its lower 
%                           and upper bounds
% 
%  - inputTrajectoryScaled: the input trajectory scaled between its lower 
%                           and upper bounds
%   
%  - transformation: the transfomation matrix for scaling (diagonal)
% 
%  - offset: the offset vector for scaling
%
%
% Example
% [x_sc, u_sc, t, c] = dataProcessing.scaleData(x, ui, lbx, ubx, lbu, ubu)
%
%
% All parameters are optional and default to zero (bounds) / empty (data).
%
% If no input trajectory is provided, then the lower and upper bounds of
% the input are ignored and only the states are scaled. If an input
% trajectory is provided, but no lower and upper bounds for the input; or 
% if the provided lower and upper bounds are equal to each other, then
% the input trajectory is not scaled and returned unchanged, the
% transformation matrix will be identity and the offset vector is zero for 
% the corresponding values. Values for states are always listed first, 
% followed by values corresponding to inputs. 
%
% Likewise for missing state trajectories or lower/upper bounds. 

% 22.02.2022 leona.schnelle@haw-hamburg.de
% February 2025 adapted by Sarah Casura 

% Define arguments and their defaults
arguments
         stateTrajectory = zeros(0)
         inputTrajectory = zeros(0)
         stateLowerBound = 0 
         stateUpperBound = 0
         inputLowerBound = 0
         inputUpperBound = 0
end

% calculate the scale factor and offset needed to scale data between 
% bounds (I moved this into a helper function to avoid code duplication). 
[stateScaleFactor, stateOffset] = dataProcessing.calculateScaleOffset( ...
    stateTrajectory, stateLowerBound, stateUpperBound);
[inputScaleFactor, inputOffset] = dataProcessing.calculateScaleOffset( ...
    inputTrajectory, inputLowerBound, inputUpperBound);

% calculate the transformation matrix and offset vector 
transformation = diag([stateScaleFactor, inputScaleFactor]);
offset = [stateOffset, inputOffset];

% calculate the scaled trajectories
stateTrajectoryScaled = (1./stateScaleFactor) .* ...
    (stateTrajectory - stateOffset); 
inputTrajectoryScaled = (1./inputScaleFactor) .* ...
    (inputTrajectory - inputOffset);
