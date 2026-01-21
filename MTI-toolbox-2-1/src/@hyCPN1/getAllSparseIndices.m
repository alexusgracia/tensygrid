function H = getAllSparseIndices(H)
% Torben Warnecke - 11/06/2024

H.sparseIndices.state = H.getSparseIndices(H.F.state);
H.sparseIndices.stateDerivative = H.getSparseIndices(H.F.stateDerivative);
H.sparseIndices.input = H.getSparseIndices(H.F.input);
H.sparseIndices.algebraic = H.getSparseIndices(H.F.algebraic);
H.sparseIndices.boolean = H.getSparseIndices(H.F.boolean);

H.sparseIndices.equality = H.getSparseIndices(H.phi.equality);
H.sparseIndices.inequality = H.getSparseIndices(H.phi.inequality);
end

