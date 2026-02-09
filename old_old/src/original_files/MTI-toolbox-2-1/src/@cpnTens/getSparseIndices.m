function [nZRows, nZColumns, nZSignVector] = getSparseIndices(sys)
    % Get Positions and dense vector of non zero factors  
    [sys.RowInd, sys.ColInd, sys.DataVec] = find(sys.U);
    
    % Detemine sign of nnz factors
    sys.SignVec = -sign(sys.DataVec);
    
    % Legacy return values
    nZRows =  sys.RowInd;
    nZColumns = sys.ColInd ;
    nZSignVector = sys.SignVec;
    
    % Ensure efficient memory layout of performance-relevant data
    sys.DataVec = sys.DataVec(:);
    sys.RowInd = sys.RowInd(:);          
end
