function sys = connectByIndices(sys, stateConnections, algebraicConnections, booleanConnections, similarInputs)
%CONNECTBYINDICES

% Torben Warnecke - 11/06/2024

% check for doubled input indices
inputIdx = [];
if ~isempty(stateConnections)
    inputIdx = [inputIdx; stateConnections(:,1)];
end
if ~isempty(algebraicConnections)
    inputIdx = [inputIdx; algebraicConnections(:,1)];
end
if ~isempty(booleanConnections)
    inputIdx = [inputIdx; booleanConnections(:,1)];
end

if length(unique(inputIdx)) ~= length(inputIdx)
    error('Input can only be connected once! Therefore input indexes can only occur once in stateConnections, algebraicConnections and booleanConnections.')
end

if ~isempty(similarInputs)
    similarInputIdx = reshape(similarInputs,[],1);
    similarInputIdx = similarInputIdx(similarInputIdx~=0);
    if length(unique(similarInputIdx)) ~= length(similarInputIdx)
        error('Input indexes can only accur once in similarInputs.')
    end

    similarInputIdx = reshape(similarInputs(:,2:end),[],1);
    similarInputIdx = similarInputIdx(similarInputIdx~=0);
    inputIdx = [inputIdx; similarInputIdx];
    if length(unique(inputIdx)) ~= length(inputIdx)
        error('Input cannot be similar to another input (and be reduced) when it also should be connected to another variable! Therefore only the first column of similarInputs can have input indices that also accur in stateConnections, algebraicConnections and booleanConnections.')
    end
end

% similar Inputs
if ~isempty(similarInputs)
    for k = 1:size(similarInputs,1)
        idx1 = similarInputs(k,1);
        idx2 = similarInputs(k,2:end);
        idx2(idx2==0) = [];
        sys.H.F.input.t(idx1,:) = sys.H.F.input.t(idx1,:) + sum(sys.H.F.input.t(idx2,:),1);
        sys.H.F.input.f(idx1,:) = sys.H.F.input.f(idx1,:) + sum(sys.H.F.input.f(idx2,:),1);
        sys.H.F.input.c(idx1,:) = sys.H.F.input.c(idx1,:) + sum(sys.H.F.input.c(idx2,:),1);
        sys.H.F.input.t(idx2,:) = 0;
        sys.H.F.input.f(idx2,:) = 0;
        sys.H.F.input.c(idx2,:) = 0;
    end
end

% state-input connections
if ~isempty(stateConnections)
    for k = 1:size(stateConnections,1)
        if sum(stateConnections(k,:)~=0) == 2
            % one variable connected to input:
            %   replace input variable with state variable
            idxIn = sign(stateConnections(k, 1)) *stateConnections(k, 1);
            idxOut = sign(stateConnections(k, 2)) *stateConnections(k, 2);

            sys.H.F.state.c(idxOut, :) = sign(stateConnections(k, 2)) * sys.H.F.state.c(idxOut, :) + sign(stateConnections(k, 1)) * sys.H.F.input.c(idxIn, :);
            sys.H.F.state.t(idxOut, :) = sign(stateConnections(k, 2)) * sys.H.F.state.t(idxOut, :) + sign(stateConnections(k, 1)) * sys.H.F.input.t(idxIn, :);
            sys.H.F.state.f(idxOut, :) = sign(stateConnections(k, 2)) * sys.H.F.state.f(idxOut, :) + sign(stateConnections(k, 1)) * sys.H.F.input.f(idxIn, :);

            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

        else
            % multiple variables connected to input:
            %   replace input variable by algebraic variable and add
            %   connection equation
            idxIn = stateConnections(k, 1);
            connections = stateConnections(k, 2:end);
            connections = connections(connections~=0);
            idxOutSigns = sign(connections);
            idxOut = abs(connections);

            sys.H.F.algebraic.c = [sys.H.F.algebraic.c; sys.H.F.input.c(idxIn,:)];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t; sys.H.F.input.t(idxIn,:)];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f; sys.H.F.input.f(idxIn,:)];
            
            sz = size(sys.H.F.algebraic.c);
            sys.H.F.algebraic.c = [sys.H.F.algebraic.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t(sz(1), sz(2)+1) = 1;

            sz = size(sys.H.F.stateDerivative.c);
            sys.H.F.stateDerivative.c = [sys.H.F.stateDerivative.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.t = [sys.H.F.stateDerivative.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.f = [sys.H.F.stateDerivative.f, zeros([sz(1), 1+length(idxOut)])];

            sys.H.F.state.c = [sys.H.F.state.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.t = [sys.H.F.state.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.f = [sys.H.F.state.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.t(idxOut, sz(2)+2:end) = diag(ones([length(idxOut),1]));
            
            sz = size(sys.H.F.boolean.c);
            sys.H.F.boolean.c = [sys.H.F.boolean.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.t = [sys.H.F.boolean.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.f = [sys.H.F.boolean.f, zeros([sz(1), 1+length(idxOut)])];
            
            sz = size(sys.H.F.input.c);
            sys.H.F.input.c = [sys.H.F.input.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.t = [sys.H.F.input.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.f = [sys.H.F.input.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

            sz = size(sys.H.phi.equality.c);
            sys.H.phi.equality.c = [sys.H.phi.equality.c, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];
            sys.H.phi.equality.t = [sys.H.phi.equality.t, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)]), 1, -idxOutSigns];
            sys.H.phi.equality.f = [sys.H.phi.equality.f, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];

            sz = size(sys.H.phi.inequality.c);
            sys.H.phi.inequality.c = [sys.H.phi.inequality.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.t = [sys.H.phi.inequality.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.f = [sys.H.phi.inequality.f, zeros([sz(1), 1+length(idxOut)])];

            if ~isempty(sys.inputName)
                sys.algebraicName = [sys.algebraicName, sys.inputName(idxIn)];
            end
            if ~isempty(sys.inputUnit)
                sys.algebraicUnit = [sys.algebraicUnit, sys.inputUnit(idxIn)];
            end
        end
    end
end

% algebraic-input connections
if ~isempty(algebraicConnections)
    for k = 1:size(algebraicConnections,1)
        if sum(algebraicConnections(k,:)~=0) == 2
            % one variable connected to input:
            %   replace input variable with state variable
            idxIn = sign(algebraicConnections(k, 1))*algebraicConnections(k, 1);
            idxOut = sign(algebraicConnections(k, 2))*algebraicConnections(k, 2);

            sys.H.F.algebraic.c(idxOut, :) = sign(algebraicConnections(k, 2))*sys.H.F.algebraic.c(idxOut, :) + sign(algebraicConnections(k, 1))*sys.H.F.input.c(idxIn, :);
            sys.H.F.algebraic.t(idxOut, :) = sign(algebraicConnections(k, 2))*sys.H.F.algebraic.t(idxOut, :) + sign(algebraicConnections(k, 1))*sys.H.F.input.t(idxIn, :);
            sys.H.F.algebraic.f(idxOut, :) = sign(algebraicConnections(k, 2))*sys.H.F.algebraic.f(idxOut, :) + sign(algebraicConnections(k, 1))*sys.H.F.input.f(idxIn, :);

            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

        else
            % multiple variables connected to input:
            %   replace input variable by algebraic variable and add
            %   connection equation
            idxIn = algebraicConnections(k, 1);
            connections = algebraicConnections(k, 2:end);
            connections = connections(connections~=0);
            idxOutSigns = sign(connections);
            idxOut = abs(connections);

            sys.H.F.algebraic.c = [sys.H.F.algebraic.c; sys.H.F.input.c(idxIn,:)];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t; sys.H.F.input.t(idxIn,:)];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f; sys.H.F.input.f(idxIn,:)];

            sz = size(sys.H.F.algebraic.c);

            sys.H.F.algebraic.c = [sys.H.F.algebraic.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t(sz(1), sz(2)+1) = 1;
            sys.H.F.algebraic.t(idxOut, sz(2)+2:end) = diag(ones([length(idxOut),1]));

            sz = size(sys.H.F.stateDerivative.c);
            sys.H.F.stateDerivative.c = [sys.H.F.stateDerivative.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.t = [sys.H.F.stateDerivative.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.f = [sys.H.F.stateDerivative.f, zeros([sz(1), 1+length(idxOut)])];

            sys.H.F.state.c = [sys.H.F.state.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.t = [sys.H.F.state.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.f = [sys.H.F.state.f, zeros([sz(1), 1+length(idxOut)])];

            sz = size(sys.H.F.boolean.c);
            sys.H.F.boolean.c = [sys.H.F.boolean.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.t = [sys.H.F.boolean.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.f = [sys.H.F.boolean.f, zeros([sz(1), 1+length(idxOut)])];
    
            sz = size(sys.H.F.input.c);
            sys.H.F.input.c = [sys.H.F.input.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.t = [sys.H.F.input.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.f = [sys.H.F.input.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

            sz = size(sys.H.phi.equality.c);
            sys.H.phi.equality.c = [sys.H.phi.equality.c, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];
            sys.H.phi.equality.t = [sys.H.phi.equality.t, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)]), 1, -idxOutSigns];
            sys.H.phi.equality.f = [sys.H.phi.equality.f, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];

            sz = size(sys.H.phi.inequality.c);
            sys.H.phi.inequality.c = [sys.H.phi.inequality.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.t = [sys.H.phi.inequality.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.f = [sys.H.phi.inequality.f, zeros([sz(1), 1+length(idxOut)])];

            if ~isempty(sys.inputName)
                sys.algebraicName = [sys.algebraicName, sys.inputName(idxIn)];
            end
            if ~isempty(sys.inputUnit)
                sys.algebraicUnit = [sys.algebraicUnit, sys.inputUnit(idxIn)];
            end
        end
    end
end

% boolean-input connections
if ~isempty(booleanConnections)
    for k = 1:size(booleanConnections,1)
        if sum(booleanConnections(k,:)~=0) == 2
            % one variable connected to input:
            %   replace input variable with state variable
            idxIn = sign(booleanConnections(k, 1))*booleanConnections(k, 1);
            idxOut = sign(booleanConnections(k, 2))*booleanConnections(k, 2);

            sys.H.F.boolean.c(idxOut, :) = sign(booleanConnections(k, 2))*sys.H.F.boolean.c(idxOut, :) + sign(booleanConnections(k, 1))*sys.H.F.input.c(idxIn, :);
            sys.H.F.boolean.t(idxOut, :) = sign(booleanConnections(k, 2))*sys.H.F.boolean.t(idxOut, :) + sign(booleanConnections(k, 1))*sys.H.F.input.t(idxIn, :);
            sys.H.F.boolean.f(idxOut, :) = sign(booleanConnections(k, 2))*sys.H.F.boolean.f(idxOut, :) + sign(booleanConnections(k, 1))*sys.H.F.input.f(idxIn, :);

            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

        else
            % multiple variables connected to input:
            %   replace input variable by algebraic variable and add
            %   connection equation
            idxIn = booleanConnections(k, 1);
            connections = booleanConnections(k, 2:end);
            connections = connections(connections~=0);
            idxOutSigns = sign(connections);
            idxOut = abs(connections);

            sys.H.F.algebraic.c = [sys.H.F.algebraic.c; sys.H.F.input.c(idxIn,:)];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t; sys.H.F.input.t(idxIn,:)];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f; sys.H.F.input.f(idxIn,:)];

            sz = size(sys.H.F.algebraic.c);

            sys.H.F.algebraic.c = [sys.H.F.algebraic.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t = [sys.H.F.algebraic.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.f = [sys.H.F.algebraic.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.algebraic.t(sz(1), sz(2)+1) = 1;

            sz = size(sys.H.F.stateDerivative.c);
            sys.H.F.stateDerivative.c = [sys.H.F.stateDerivative.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.t = [sys.H.F.stateDerivative.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.stateDerivative.f = [sys.H.F.stateDerivative.f, zeros([sz(1), 1+length(idxOut)])];

            sys.H.F.state.c = [sys.H.F.state.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.t = [sys.H.F.state.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.state.f = [sys.H.F.state.f, zeros([sz(1), 1+length(idxOut)])];

            sz = size(sys.H.F.boolean.c);
            sys.H.F.boolean.c = [sys.H.F.boolean.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.t = [sys.H.F.boolean.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.f = [sys.H.F.boolean.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.boolean.t(idxOut, sz(2)+2:end) = diag(ones([length(idxOut),1]));

            sz = size(sys.H.F.input.c);
            sys.H.F.input.c = [sys.H.F.input.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.t = [sys.H.F.input.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.f = [sys.H.F.input.f, zeros([sz(1), 1+length(idxOut)])];
            sys.H.F.input.c(idxIn, :) = 0;
            sys.H.F.input.t(idxIn, :) = 0;
            sys.H.F.input.f(idxIn, :) = 0;

            sz = size(sys.H.phi.equality.c);
            sys.H.phi.equality.c = [sys.H.phi.equality.c, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];
            sys.H.phi.equality.t = [sys.H.phi.equality.t, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)]), 1, -idxOutSigns];
            sys.H.phi.equality.f = [sys.H.phi.equality.f, zeros([sz(1), 1+length(idxOut)]);...
                zeros([1, sz(2)+1+length(idxOut)])];

            sz = size(sys.H.phi.inequality.c);
            sys.H.phi.inequality.c = [sys.H.phi.inequality.c, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.t = [sys.H.phi.inequality.t, zeros([sz(1), 1+length(idxOut)])];
            sys.H.phi.inequality.f = [sys.H.phi.inequality.f, zeros([sz(1), 1+length(idxOut)])];

            if ~isempty(sys.inputName)
                sys.algebraicName = [sys.algebraicName, sys.inputName(idxIn)];
            end
            if ~isempty(sys.inputUnit)
                sys.algebraicUnit = [sys.algebraicUnit, sys.inputUnit(idxIn)];
            end
        end
    end
end

% delete empty input matrices
sz = size(sys.H.F.input.c);
[rowC, ~] = find(sys.H.F.input.c);
[rowT, ~] = find(sys.H.F.input.t);
[rowF, ~] = find(sys.H.F.input.f);
emptyIdx = ~any([rowC; rowT; rowF] == 1:sz(1), 1);
sys.H.F.input.c(emptyIdx, :) = [];
sys.H.F.input.t(emptyIdx, :) = [];
sys.H.F.input.f(emptyIdx, :) = [];


sys.m = size(sys.H.F.input.c,1);
sys.p = size(sys.H.F.algebraic.c,1);

if ~isempty(sys.inputName)
    sys.inputName(emptyIdx) = [];
end
if ~isempty(sys.inputUnit)
    sys.inputUnit(emptyIdx) = [];
end

for k = 1:length(emptyIdx)
    sys.manipulatedVariable(find(sys.manipulatedVariable==emptyIdx(k))) = [];
    sys.measuredDisturbance(find(sys.measuredDisturbance==emptyIdx(k))) = [];
    sys.unmeasuredDisturbance(find(sys.unmeasuredDisturbance==emptyIdx(k))) = [];
end


end

