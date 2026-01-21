function [scaleFactor, offset] = ...
    calculateScaleOffset(trajectory, lowerBound, upperBound)
%
% <calculateScaleOffset> - Helper function to find scale factor and offset 
% given data trajectory (to avoid code duplication in scaleData.m). 
% 
% % Input parameters: 
%
%  - trajectory: data (typically state or input trajectory)
%                This can be a matrix containing several trajectories for
%                several states/inputs (one trajectory per column). 
%
%  - lowerBound: lower bound(s) to scale the data to
%                This can be a scalar (in which case it is applied to all
%                trajectories), or a vector with the number of elements
%                corresponding to the number of columns in trajectory. 
%
%  - upperBound: upper bound(s) to scale the data to
%                This can be a scalar (in which case it is applied to all
%                trajectories), or a vector with the number of elements
%                corresponding to the number of columns in trajectory. 
%
% Output parameters: 
%
%  - scaleFactor: the scale factor(s) needed to scale data between bounds
%                 (defined such that data needs to be divided by it). 
% 
%  - offset: the offset(s) needed to scale data between bounds
%            (defined such that it needs to be subtracted from data AFTER
%             dividing by the scale factor). 
%
%
% Example
% [scaleFactor, offset] = dataProcessing.calculateScaleOffset( ...
%    trajectory, lowerBound, upperBound);
%
%
% All parameters are optional and default to zero (bounds) / empty (data).
%
% If no trajectory is provided, then scaleFactor and offset will be empty. 
% If a trajectory is provided but no bounds, then the scaleFactor will be 
% 1 and the offset will be 0 (corresponding to no scaling). The same
% applies if the lower and upper bounds equal each other. 

% Sarah Casura, February 2025 (adapted from scaleData.m)

% Define arguments and their defaults
arguments
         trajectory = zeros(0);
         lowerBound = 0; 
         upperBound = 0;
         % turns out the MATLAB name-value syntax is very limited in its
         % use. So it's not actually possible to call the function like in
         % any other language: calculateScaleOffset(trajectory = x, ...
         % lowerBound = 1, upperBound = 2) errors with "too many input
         % arguments". calculateScaleOffset(x, 1, 2) works fine. Hence,
         % calculateScaleOffset(x) or calculateScaleOffset() works as
         % intended (it just takes the defaults), but there is no way to
         % easily specify only bounds and no trajectory. Workarounds would
         % be to use inputParser or to use a structure, e.g. having only
         % one input called options; and then specifying arguments as
         % options.trajectory, options.lowerBound, options.upperBound. Then
         % it works with the name-value syntax. But that makes the function
         % much less readable and is not worth the effort given that this
         % is very much a corner case in the first place. 
end

% calculate scale factor and offset 
if ~isempty(trajectory)

    % make sure all bounds are the correct length 
    [~, numberOfBoundsExpected] = size(trajectory);
    if length(lowerBound) ~= numberOfBoundsExpected 
        if(length(lowerBound) ~= 1)
            warning("Lower bound must be scalar or of the same " + ...
                "length as the number of columns in trajectory.\n" + ... 
                "Applying the first lower bound to all trajectories.")
            % print a warning if it is the wrong length; and use the first
            % value of lowerBound only. Then make scalar into vector. 
        end
        lowerBound = zeros(1, numberOfBoundsExpected) + lowerBound(1);
    end
    if length(upperBound) ~= numberOfBoundsExpected 
        if(length(upperBound) ~= 1)
            warning("Upper bound must be scalar or of the same " + ...
                "length as the number of columns in trajectory.\n" + ... 
                "Applying the first upper bound to all trajectories.")
            % print a warning if it is the wrong length; and use the first
            % value of lowerBound only. Then make scalar into vector. 
        end
        upperBound = zeros(1, numberOfBoundsExpected) + upperBound(1);
    end
    flippedBounds = lowerBound > upperBound; % check for flipped bounds
    if any(flippedBounds) % flip the bounds
        warning("Upper bound is smaller than lower bound. Flipping them.")
        tempBound = lowerBound(flippedBounds);
        lowerBound(flippedBounds) = upperBound(flippedBounds);
        upperBound(flippedBounds) = tempBound;
    end 
    % scale the data to be btw. lower and upper bound
    scaleFactor = (max(trajectory) - min(trajectory)) ...
        ./ (upperBound - lowerBound);
    % in case the trajectory contains all the same values, we get
    % zero scale factor, which we change to 1, i.e. no scaling. 
    scaleFactor(scaleFactor == 0) = 1;
    offset = min(trajectory) - lowerBound .* scaleFactor;
    % in case the lower and upper bounds equal each other, we get infinite
    % scale factor and offset, which we change to 1 and 0 (no scaling)
    scaleFactor(upperBound == lowerBound) = 1; 
    offset(upperBound == lowerBound) = 0; 
else % trajectory was empty, so leave scale factor and offset empty
    scaleFactor = zeros(0);
    offset = zeros(0); 
end

