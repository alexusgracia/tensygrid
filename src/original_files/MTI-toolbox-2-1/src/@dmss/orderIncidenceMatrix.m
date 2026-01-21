function [Isorted, rowSort, colSort] = orderIncidenceMatrix(sys,I)
% Torben Warnecke - 11/06/2024

[n,m] = size(I);
if n~=m
    warning('The system is not a square system! The number of unknown variables do not equal the number of equations and constraints.')
end

% SORT EQUATIONS IN ORDER OF THE NUMBER OF VARIABLES IN THEM
% (DUMALGE-MENDELSOHN DECOMPOSITION)
[p,q,r,s,cc,rr] = dmperm(I);

Isorted = I(flip(p), flip(q));
rowSort = flip(p)';
colSort = flip(q)';

end

