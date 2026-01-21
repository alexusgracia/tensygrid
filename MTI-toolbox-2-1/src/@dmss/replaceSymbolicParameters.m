function sys = replaceSymbolicParameters(sys, old, new)
%REPLACESYMBOLICPARAMETERS Replace symbolic parameter values of the model tensor by new values
% (numeric or symbolic).
%
%   input Paramter:
%   - sys: dmss model
%   - old: old symbolic parameter values as symbolic array, which will be
%   replaced
%   - new: new parameter values either symbolic or numeric array with the
%   same length as old
%
% For detailed documentation see <a href="matlab:open((which('replaceSymbolicParametersDoc.html')))">here</a>

% Torben Warnecke - 11/06/2024
arguments
    sys dmss
    old (1,:) sym
    new (1,:)
end
sys.H = replaceSymbolicValue(sys.H, old, new);

try 
    sys.H = sys.H.convert2doubles;
catch 
end

end

