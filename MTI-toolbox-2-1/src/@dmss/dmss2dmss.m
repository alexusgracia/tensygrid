function sys = dmss2dmss(sys, stateOffset, stateSlope, inputOffset, inputSlope, algebraicOffset, algebraicSlope)
%DMSS2DMSS coordinate transformation of variables for dmss
%models defined by offsets and slopes (x = A*xT + B)
%
%   input Paramter:
%   - sys: dmss-object which coordinates will be transformed.
%   - stateOffset: Offset vector for linear state coordinate
%   transformations.
%   - stateSlope: Slope vector (or diagonal Transformation matrix) for linear state coordinate
%   transformations. Nondiagonal transformation matrices not implemented
%   yet.
%   - inputOffset: Offset vector for linear input coordinate
%   transformations.
%   - inputSlope: Slope vector (or diagonal Transformation matrix) for linear input coordinate
%   transformations. Nondiagonal transformation matrices not implemented
%   yet.
%   - algebraicOffset: Offset vector for linear algebraic variable coordinate
%   transformations.
%   - algebraicSlope: Slope vector (or diagonal Transformation matrix) for linear algebraic variable coordinate
%   transformations. Nondiagonal transformation matrices not implemented
%   yet.
%
%   output Parameter:
%   - sys: transformed system as dmss-object

% Torben Warnecke - 11/06/2024

arguments
    sys
    stateOffset         (:, 1) double
    stateSlope          (:, :) double
    inputOffset         (:, 1) double
    inputSlope          (:, :) double
    algebraicOffset     (:, 1) double
    algebraicSlope      (:, :) double
end

if (~isempty(stateOffset))
    if (length(stateOffset) ~= sys.n)
        error('Length of the stateOffset vector must equal the number n of states of the model.')
    end
end
if (~isempty(stateSlope))
    [n,m] = size(stateSlope);
    if (n == sys.n) && (m == sys.n)
        if sum(any(stateSlope)) < n
            error('transformed system does not represent a square system. Transformation wont be computed. Make sure that the transformation matrices do not have any all zeros rows or columns.')
        elseif ~isdiag(stateSlope)
            error('Functionality for non diagonal transformation matrices not implemented yet, since the system will in general result in a polynomial system.')
        else
            stateSlope = diag(stateSlope);
        end
    elseif (n == sys.n) && (m == 1)
    elseif (n == 1) && (m == sys.n)
        stateSlope = stateSlope';
    else
        error('stateSlope must be a vector or square matrix with the size equal to the number n of states of the model.')
    end
end

if (~isempty(inputOffset))
    if (length(inputOffset) ~= sys.m)
        error('Length of the inputOffset vector must equal the number m of inputs of the model.')
    end
end
if (~isempty(inputSlope))
    [n,m] = size(inputSlope);
    if (n == sys.m) && (m == sys.m)
        if sum(any(inputSlope)) < n
            error('transformed system does not represent a square system. Transformation wont be computed. Make sure that the transformation matrices do not have any all zeros rows or columns.')
        elseif ~isdiag(inputSlope)
            error('Functionality for non diagonal transformation matrices not implemented yet, since the system will in general result in a polynomial system.')
        else
            inputSlope = diag(inputSlope);
        end
    elseif (n == sys.m) && (m == 1)
    elseif (n == 1) && (m == sys.m)
        inputSlope = inputSlope';
    else
        error('inputSlope must be a vector or square matrix with the size equal to the number m of inputs of the model.')
    end
end

if (~isempty(algebraicOffset))
    if (length(algebraicOffset) ~= sys.p)
        error('Length of the algebraicOffset vector must equal the number p of algebraic variables of the model.')
    end
end
if (~isempty(algebraicSlope))
    [n,m] = size(algebraicSlope);
    if (n == sys.p) && (m == sys.p)
        if sum(any(algebraicSlope)) < n
            error('transformed system does not represent a square system. Transformation wont be computed. Make sure that the transformation matrices do not have any all zeros rows or columns.')
        elseif ~isdiag(algebraicSlope)
            error('Functionality for non diagonal transformation matrices not implemented yet, since the system will in general result in a polynomial system.')
        else
            algebraicSlope = diag(algebraicSlope);
        end
    elseif (n == sys.p) && (m == 1)
    elseif (n == 1) && (m == sys.p)
        algebraicSlope = algebraicSlope';
    else
        error('algebraicSlope must be a vector or square matrix with the size equal to the number m of algerbaic variables of the model.')
    end
end

% Convert all symbolic values to double values (currently needed)
sys.H = sys.H.convert2doubles();

% principle:
%   [1-|F| + F*x] = [1-|F| + F*(a * xtilde + b)]
%   = (1-|F|+F*(a+b)) * [1-|F*a/(1-|F|+F*(a+b))| + F*a/(1-|F|+F*(a+b)) * xtilde]

phiEq = sys.H.phi.equality.c + sys.H.phi.equality.t - sys.H.phi.equality.f;
phiIneq = sys.H.phi.inequality.c + sys.H.phi.inequality.t - sys.H.phi.inequality.f;

    function [HFnew, phiEq, phiIneq] = variableTransformation(HFold, phiEq, phiIneq, slope, offset)
        if isempty(offset)
            offset = zeros([length(slope),1]);
        end
        if isempty(slope)
            slope = ones([length(slope),1]);
        end
        Hphinewc = abs(1-abs(HFold.c) + HFold.c .* offset) + abs(HFold.c .* slope);
        Hphinewt = (abs(slope)+abs(offset));
        Hphinewf = 2*(abs(1/2-1/2*offset) + abs(-1/2*slope));

        sgnc = 2*((1-abs(HFold.c) + HFold.c .* offset)>=0) -1;
        sgnt = 2*((offset)>=0) -1;
        sgnf = 2*((1/2-1/2*offset)>=0) -1;

        HFnewc = HFold.c .* slope ./ Hphinewc .* sgnc;
        HFnewt = HFold.t .* slope ./ Hphinewt .* sgnt;
        HFnewf = -HFold.f.* slope ./ Hphinewf .* sgnf;

        HFnew = HFnewc + HFnewt + HFnewf;

        Hphinew = (HFold.c~=0) .* Hphinewc.*sgnc + HFold.t .* Hphinewt .* sgnt + HFold.f .* Hphinewf .* sgnf;
        Hphinew(Hphinew==0) = 1;

        phiEq = phiEq .* prod(Hphinew, 1);
        if ~isempty(phiIneq)
            phiIneq = phiIneq .* prod(Hphinew, 1);
        end
    end

%% state transformation
if (~isempty(stateOffset) || ~isempty(stateSlope))
    [HFx, phiEq, phiIneq] = variableTransformation(sys.H.F.state, phiEq, phiIneq, stateSlope, stateOffset);
end

%% state derivative transformation
if (~isempty(stateOffset) || ~isempty(stateSlope))
    if sys.ts > 0
        [HFxp, phiEq, phiIneq] = variableTransformation(sys.H.F.stateDerivative, phiEq, phiIneq, stateSlope, stateOffset);
    elseif sys.ts == 0
        [HFxp, phiEq, phiIneq] = variableTransformation(sys.H.F.stateDerivative, phiEq, phiIneq, stateSlope, []);
    end
end

%% algebraic variable transformation
if (~isempty(algebraicOffset) || ~isempty(algebraicSlope))
    [HFy, phiEq, phiIneq] = variableTransformation(sys.H.F.algebraic, phiEq, phiIneq, algebraicSlope, algebraicOffset);
end

%% input transformation
if (~isempty(inputOffset) || ~isempty(inputSlope))
    [HFu, phiEq, phiIneq] = variableTransformation(sys.H.F.input, phiEq, phiIneq, inputSlope, inputOffset);
end

%% 
if (~isempty(stateOffset) || ~isempty(stateSlope))
    sys.H.F.state = HFx;
    sys.H.F.stateDerivative = HFxp;
end
if (~isempty(algebraicOffset) || ~isempty(algebraicSlope))
    sys.H.F.algebraic = HFy;
end
if (~isempty(inputOffset) || ~isempty(inputSlope))
    sys.H.F.input = HFu;
end

sys.H.phi.equality = phiEq;
sys.H.phi.inequality = phiIneq;

end

