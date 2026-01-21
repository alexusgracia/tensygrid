function [sys, stateOffset, stateSlope, inputOffset, inputSlope, algebraicOffset, algebraicSlope] = normalize(sys,stateCurrentRange, stateTargetRange, inputCurrentRange, inputTargetRange, algebraicCurrentRange, algebraicTargetRange)
%NORMALIZE Coordinate transformation of variables for dmss
%models defined by old and new ranges
%
%   input Paramter:
%   - sys: dmss-object which coordinates will be transformed.
%   - stateCurrentRange: Range of the states that should be transformed
%   into the target Range.
%   - stateTargetRange: Target range for the state transformation (e.g. [0,
%   1]).
%   - inputCurrentRange: Range of the inputs that should be transformed
%   into the target Range.
%   - inputTargetRange: Target range for the input transformation.
%   - algebraicCurrentRange: Range of the algebraic variables that should 
%   be transformed into the target Range.
%   - algebraicTargetRange: Target range for the algebraic variable
%   transformation.
%
%   output Parameter:
%   - sys: transformed system as dmss-object
%   - stateOffset: Offset vector for linear state coordinate
%   transformations.
%   - stateSlope: Slope vector for linear state coordinate transformations.
%   - inputOffset: Offset vector for linear input coordinate
%   transformations.
%   - inputSlope: Slope vector for linear input coordinate transformations. 
%   - algebraicOffset: Offset vector for linear algebraic variable coordinate
%   transformations.
%   - algebraicSlope: Slope vector for linear algebraic variable coordinate
%   transformations.
%
% For detailed documentation see <a href="matlab:open((which('normalizeDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024

arguments
    sys
    stateCurrentRange         (:, 2) double
    stateTargetRange          (:, 2) double
    inputCurrentRange         (:, 2) double
    inputTargetRange          (:, 2) double
    algebraicCurrentRange     (:, 2) double
    algebraicTargetRange      (:, 2) double
end

if isempty(stateCurrentRange) && isempty(stateTargetRange)
    stateOffset = [];
    stateSlope = [];
elseif (sum(size(stateCurrentRange)~=[sys.n,2])>0)||(sum(size(stateTargetRange)~=[sys.n,2])>0)
    error('Size of the current and Target state Range matrices must be [n,2], with the minimal values as first columen and maximal values as second column. Or both can be empty, then no transformation will be performed.')
else
    xmin = stateCurrentRange(:,1);
    xmax = stateCurrentRange(:,2);
    xTmin = stateTargetRange(:,1);
    xTmax = stateTargetRange(:,2);

    stateSlope = (xmax-xmin)'./(xTmax-xTmin)';
    stateOffset = xmin' - (stateSlope.*xTmin');
end

if isempty(inputCurrentRange) && isempty(inputTargetRange)
    inputOffset = [];
    inputSlope = [];
elseif (sum(size(inputCurrentRange)~=[sys.m,2])>0)||(sum(size(inputTargetRange)~=[sys.m,2])>0)
    error('Size of the current and Target input Range matrices must be [m,2], with the minimal values as first columen and maximal values as second column. Or both can be empty, then no transformation will be performed.')
else
    xmin = inputCurrentRange(:,1);
    xmax = inputCurrentRange(:,2);
    xTmin = inputTargetRange(:,1);
    xTmax = inputTargetRange(:,2);

    inputSlope = (xmax-xmin)'./(xTmax-xTmin)';
    inputOffset = xmin' - (inputSlope.*xTmin');
end

if isempty(algebraicCurrentRange) && isempty(algebraicTargetRange)
    algebraicOffset = [];
    algebraicSlope = [];
elseif (sum(size(algebraicCurrentRange)~=[sys.p,2])>0)||(sum(size(algebraicTargetRange)~=[sys.p,2])>0)
    error('Size of the current and Target algebraic variable Range matrices must be [p,2], with the minimal values as first columen and maximal values as second column. Or both can be empty, then no transformation will be performed.')
else
    xmin = algebraicCurrentRange(:,1);
    xmax = algebraicCurrentRange(:,2);
    xTmin = algebraicTargetRange(:,1);
    xTmax = algebraicTargetRange(:,2);

    algebraicSlope = (xmax-xmin)'./(xTmax-xTmin)';
    algebraicOffset = xmin' - (algebraicSlope.*xTmin');
end

sys = dmss2dmss(sys, stateOffset, stateSlope, inputOffset, inputSlope, algebraicOffset, algebraicSlope);

end

