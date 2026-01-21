import numpy as np
import scipy as sp
import time

import memory_profiler as mp

@mp.profile
def main():
    # %% Start file for TenSyGrid Hackathon 16.12.2025
    # Linearization for iMTI models in CPN representation (basic, functional)

    # Configuration
    check = 0  # 0 = sparse/scalable/fast | 1 = full (only for checks)
    test = False
    test_zeros = False
    debug = True
    # %% Define problem


    # iMTI model in CPN representation
    if test:
        n = 2               # number of states
        m = 1               # number of inputs
        p = 0               # number of outputs
        q = n               # number of equations
        N = 2 * n + m + p   # total number of signals

        S = np.array([
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 1],
            [0, 0, 0, 1, 0, 1],
            [0, 0, 0, 0, 1, 0],
        ]) # Structure matrix
        P = np.array([
            [-1, 0, -1, 0, 1, 0],
            [0, -1, 0, -1, 0, 1],
        ]) # Parameter matrix
        r = S.shape[1] # rank (number of non-zero elements in the structure matrix)
        S= sp.sparse.csr_matrix(S)
        P = sp.sparse.csr_matrix(P)

        # Linearization point
        if not test_zeros:
            dx = np.ones((n, 1))    # state derivative 
        else:
            dx = np.zeros((n, 1))   # (CONFLICTIVE CASE) state derivative (TODO: Develop)
        x = 2*np.ones((n, 1))       # state
        u = np.ones((m, 1))         # input
        y = np.ones((p, 1))         # output

    else:
        n = 1000                    # number of states
        m = 1                       # number of inputs
        p = 0                       # number of outputs
        q = n                       # number of equations
        N = 2 * n + m + p           # total number of signals
        r = 1000                       # rank (number of non-zero elements in the structure matrix)
        S = np.random.rand(N, r)    # Structure matrix
        S = np.round(2*S-1)
        P = np.random.rand(q, r)    # Parameter matrix
        P = P*(P<0.3)
        S = sp.sparse.csr_matrix(S)
        P = sp.sparse.csr_matrix(P)

        # linearization point
        dx = np.random.rand(n, 1)   # state derivative 
        x = np.random.rand(n, 1)    # state
        u = np.random.rand(m, 1)    # input
        y = np.random.rand(p, 1)    # output

    # signal vector
    v = np.vstack([dx, x, u, y])        # signal vector

    tic = time.process_time()   # start the cputime clock
    # %% Linearization (Scalable = sparse & low rank)
    if debug: print(S.toarray())
    X = S.multiply(v) - abs(S)
    if debug: print(X.toarray())
    icrit = X == -1 # get numerical critical indices
    if icrit.nnz > 0:
        raise RuntimeError("specials to be implemented ?")

    X.data = X.data + 1.0 # add one only to the nonzero elements

    # We want, for all S != 0 and X == 0, to set X to 1 efficiently (without full matrix loops).
    # Both S and X are sparse.
    # So: find the (i,j) where S is nonzero but X is 0, and efficiently set X[i,j]=1.
    S_nonzero = S.nonzero()
    X_nz_dict = set(zip(*X.nonzero()))

    # Collect new (i,j) that are in S but *not* in X (i.e., S!=0 and X==0)
    to_add = []
    for i, j in zip(*S_nonzero):
        if (i, j) not in X_nz_dict:
            to_add.append((i, j))

    if to_add:  # if there are positions to modify
        # prepare data to add (location and value = 1.0)
        rows, cols = zip(*to_add)
        data = np.ones(len(to_add))
        shape = X.shape
        X_add = sp.sparse.coo_matrix((data, (rows, cols)), shape=shape)
        X = X + X_add
        X = X.tocsr()  # ensure CSR format for downstream code

    if debug: print(X.toarray())

    # get the indices and values of X
    rowi, coli = X.nonzero()
    val = X.data

    Y = np.zeros(r) # initialize product vector
    for c in range(r): # product over nonzero elements
        mask = coli == c # get the indices of the non-zero elements in the column c
        if np.any(mask): # if there are non-zero elements in the column c
            Y[c] = np.prod(val[mask]) # compute the product of the non-zero elements in the column c

    X.data = 1.0 / X.data # Invert only nonzero elements of X
    F = S.multiply(Y).multiply(X)

    EABC = P*(F.T) 
    E = -EABC[:, 0:n] # state matrix
    A = EABC[:, n:2*n] # system matrix
    B = EABC[:, 2*n:2*n+m] # input matrix

    # Local stability: dominant (geralized) eigenvalues
    A_dense = A.toarray()
    E_dense = E.toarray()

    #eigvals, eigvecs = sp.linalg.eig(A_dense, E_dense)  # resuelve A v = λ E v
    #print(eigvals)
    toc = time.process_time()
    print( "Time: ", toc-tic)

if __name__ == "__main__":
    mp.profile(main())