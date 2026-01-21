function T = symSingleVars(symbol,N)
% Function symSingleVars creates a vector of symbolic variables with the given name
% and dimension.
    T = sym(zeros(N,1)); 
    for k = 1:N
          syms(sprintf('%s%d(t)', symbol,k))
          TT = symfun(eval(sprintf('%s%d(t)', symbol,k)), t); % using eval instead of sym
          T(k,:) = TT;
    end
end