function Out = getSparseIndices(~,F)
% Torben Warnecke - 11/06/2024

nCols = size(F.t,2);
[tRow, tCol] = find(F.t);
[fRow, fCol] = find(F.f);
%[tfRow, tfCol, tfData] = find(F.t-F.f);
%tfOne = tfData>0;
[cRow, cCol, cData] = find(F.c);
cSigns = -sign(cData);

Out.t = [reshape(tRow,[],1), reshape(tCol,[],1)];
Out.f = [reshape(fRow,[],1), reshape(fCol,[],1)];
Out.c = [reshape(cRow,[],1), reshape(cCol,[],1), reshape(cData,[],1), reshape(cSigns,[],1)];
Out.nCols = nCols;
end

