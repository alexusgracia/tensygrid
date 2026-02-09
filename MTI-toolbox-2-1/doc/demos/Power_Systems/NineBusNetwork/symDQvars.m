function T = symDQvars(symbol,N)
    % Creates N pairs of dq-variables of the given name
    v = ['d','q'];
    T = sym(zeros(2*N,1)); 
    for k = 0:N-1
        for i=1:2
          syms(sprintf('%s%d%s(t)', symbol,k+1,v(i)))
          TT = symfun(eval(sprintf('%s%d%s(t)', symbol,k+1,v(i))), t); % use eval instead of sym
          T(i+2*k) = TT;
        end
    end
end