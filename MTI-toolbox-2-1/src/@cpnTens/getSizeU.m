function [rows, cols] = getSizeU(sys)
    [sys.colsU, sys.rowsU] = size(sys.U);
    rows = sys.rowsU;
    cols = sys.colsU;
end
