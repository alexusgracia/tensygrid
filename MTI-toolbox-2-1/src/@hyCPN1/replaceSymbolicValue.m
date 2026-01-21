function H = replaceSymbolicValue(H, old, new)
%REPLACESYMBOLICVALUES 

% Torben Warnecke - 11/06/2024
arguments
    H
    old (1,:) sym
    new (1,:)
end

H.F.stateDerivative.c = subs(H.F.stateDerivative.c, old, new);
H.F.state.c = subs(H.F.state.c, old, new);
H.F.input.c = subs(H.F.input.c, old, new);
H.F.algebraic.c = subs(H.F.algebraic.c, old, new);
H.F.boolean.c = subs(H.F.boolean.c, old, new);

H.phi.equality.c = subs(H.phi.equality.c, old, new);
H.phi.inequality.c = subs(H.phi.inequality.c, old, new);

try
    H = H.convert2doubles();
catch
    disp('Info: There are still symbolic parameters in the model tensor.')
end

end

