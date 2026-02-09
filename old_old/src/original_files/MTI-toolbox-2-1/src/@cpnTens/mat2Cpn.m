function [U,phi,rank] = mat2Cpn(M,normtype)
%MAT2CPN Direct converstion of the matrix representation of a multilinear 
%   tensor to CPN1.
%  M must be a matrix of n x 2^r dimensions, with r being an integer.
    arguments
        M
        normtype char = '1'
    end

    if(normtype ~= '1')
        error("Only normtype '1' implemented as of now.")
    end

    N = log2(size(M,2));
    if N~=round(N)
        error("Second dimension of the matrix must be 2^r, with r being an integer.")
    end

    U = zeros([N,1]);
    for k = 1:N
        Uk = U;
        Uk(k,:) = 1;
        U = [U, Uk];
    end
    phi = M;
    rank = size(M,2);
end

