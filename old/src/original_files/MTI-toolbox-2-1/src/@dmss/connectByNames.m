function sys = connectByNames(sys)
%CONNECTBYNAMES
% creates connection matrices based of the names

% Torben Warnecke - 11/06/2024

stateConnections = [];
algebraicConnections = [];
booleanConnections = [];
similarInputs = [];

% if multiple inputs have the same name, they will be converted into one input
%[~, idxRed, idxFull] = unique(sys.inputName,'stable');

sz = 1;
idx1 = 1;
nSimilar = 0;
while idx1 < length(sys.inputName)
    if ~any(similarInputs(:,2:end)==idx1,'all')
        idxInput = find(matches(sys.inputName(idx1+1:end), sys.inputName(idx1)));
        if ~isempty(idxInput)
            m = length(idxInput);
            if m>sz && nSimilar~=0
                similarInputs(:, 1+sz+1:1+m) = 0;
                sz = m;
            end
            similarInputs(end+1,1:1+m) = [idx1, idxInput+idx1];
            nSimilar = nSimilar+1;
        end
    end
    idx1 = idx1+1;
end
if ~isempty(similarInputs)
    disp('Info: Multiple Inputs have the same name. They will be converted into one similar input!')
    sys = sys.connectByIndices([], [], [], similarInputs);
end


% connecting inputs and signals (states, algebraic, boolean)
connectedNames = [];

% states
idxInput = 1;
while idxInput <= length(sys.inputName)
    idxOut = find(matches(sys.stateName, sys.inputName(idxInput)));
    if length(idxOut) > 1
        error('Multiple state variables have the same name. Use unique names.')
    elseif ~isempty(idxOut)
        stateConnections = [stateConnections; idxInput, idxOut];
        connectedNames = [connectedNames, sys.stateName(idxOut)];
    end
    idxInput = idxInput+1;
end

% algebraic
idxInput = 1;
while idxInput <= length(sys.inputName)
    idxOut = find(matches(sys.algebraicName, sys.inputName(idxInput)));
    if length(idxOut) > 1
        error('Multiple algebraic variables have the same name. Use unique names.')
    elseif ~isempty(idxOut)
        algebraicConnections = [algebraicConnections; idxInput, idxOut];
        connectedNames = [connectedNames, sys.algebraicName(idxOut)];
    end
    idxInput = idxInput+1;
end

% boolean
idxInput = 1;
while idxInput <= length(sys.inputName)
    idxOut = find(matches(sys.booleanName, sys.inputName(idxInput)));
    if length(idxOut) > 1
        error('Multiple boolean variables have the same name. Use unique names.')
    elseif ~isempty(idxOut)
        booleanConnections = [booleanConnections; idxInput, idxOut];
        connectedNames = [connectedNames, sys.booleanName(idxOut)];
    end
    idxInput = idxInput+1;
end

if isempty(connectedNames)
    warning('There are no similar signal names. No signals have been connected.')
else
    for k  = 1:length(connectedNames)
        if sum(matches(connectedNames(k+1:end), connectedNames(k))) > 0
            error('Multiple output variables (state and algebraic variables) have the same name. Use unique names for output signals.')
        end
    end
    sys = sys.connectByIndices(stateConnections, algebraicConnections, booleanConnections, []);
end

end

