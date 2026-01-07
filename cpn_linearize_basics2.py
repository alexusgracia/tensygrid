import numpy as np
import scipy as sp
import time

# %% Start file for TenSyGrid Hackathon 16.12.2025
# Linearization for iMTI models in CPN representation (basic, functional)

# Configuration
check = 0  # 0 = sparse/scalable/fast | 1 = full (only for checks)
test = False
test_zeros = False
debug = False
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
    n = 2                       # number of states
    m = 1                       # number of inputs
    p = 0                       # number of outputs
    q = n                       # number of equations
    N = 2 * n + m + p           # total number of signals
    r = 6                       # rank (number of non-zero elements in the structure matrix)
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
X = S.multiply(v) - abs(S)  # sparse matrix of all Factors - 1
if debug: print(X.toarray())

# Compute critical indices efficiently using sparse operations
# izero: indices where S.*v == abs(S) and S != 0 (i.e., where X == 0 but S is nonzero)
# Efficiency: This is more efficient than checking all positions - we only check where S is nonzero
Sv = S.multiply(v)
absS = abs(S)
# Use dense comparison for accuracy (sparse comparison can miss zero values)
# Efficiency: only convert to dense for comparison, then back to sparse
Sv_dense = Sv.toarray()
absS_dense = absS.toarray()
S_dense = S.toarray()
izero_mask = (np.abs(Sv_dense - absS_dense) < 1e-10) & (S_dense != 0)
izero = sp.sparse.csr_matrix(izero_mask.astype(float))

# imone: indices where X == -1 (critical case before adding 1)
# Efficiency: sparse boolean operation is fast, only stores nonzero positions
imone = sp.sparse.csr_matrix((X == -1).astype(float))

# ctwo: columns where more than one X == -1
# Efficiency: count occurrences per column directly from sparse indices, O(nnz) operation
imone_rows, imone_cols = imone.nonzero()
# Count how many times each column appears (i.e., how many X==-1 per column)
if len(imone_cols) > 0:
    # Use bincount to count occurrences per column, pad to size r
    col_counts = np.bincount(imone_cols, minlength=r)
    ctwo = col_counts > 1
else:
    # No X == -1 found, so no problematic columns
    ctwo = np.zeros(r, dtype=bool)

# Check for critical cases (X == -1) before proceeding
icrit = X == -1  # get numerical critical indices
if icrit.nnz > 0:
    raise RuntimeError("specials to be implemented ?")

# Add one to all nonzero elements (efficient: operates only on data array, not full matrix)
X.data = X.data + 1.0

# Correct both izero and imone cases in one operation
# Efficiency: combine corrections to avoid multiple sparse matrix operations
# X(izero)=1 corrects indices in S but not X and X(imone)=1 for one X=-1
izero_coo = izero.tocoo()
imone_coo = imone.tocoo()
X_coo = X.tocoo()

# Get all positions that need correction (union of izero and imone)
# Efficiency: use set operations for O(1) lookup instead of nested loops
izero_positions = set(zip(izero_coo.row, izero_coo.col))
imone_positions = set(zip(imone_coo.row, imone_coo.col))
all_correction_positions = izero_positions | imone_positions

if all_correction_positions:
    # Create mapping for efficient position lookup
    # Efficiency: dictionary lookup O(1) vs linear search O(n)
    X_positions = {}
    for idx, (row_idx, col_idx) in enumerate(zip(X_coo.row, X_coo.col)):
        X_positions[(row_idx, col_idx)] = idx
    
    # Lists to store new positions if needed
    new_rows = []
    new_cols = []
    new_data = []
    
    # Update existing positions and collect new ones
    for (i, j) in all_correction_positions:
        if (i, j) in X_positions:
            X_coo.data[X_positions[(i, j)]] = 1.0
        else:
            new_rows.append(i)
            new_cols.append(j)
            new_data.append(1.0)
    
    # Add new positions if any (efficient: concatenate arrays once)
    if new_rows:
        X_coo.row = np.concatenate([X_coo.row, new_rows])
        X_coo.col = np.concatenate([X_coo.col, new_cols])
        X_coo.data = np.concatenate([X_coo.data, new_data])
    
    X = X_coo.tocsr()

if debug: print(X.toarray())

# get the indices and values of X
rowi, coli = X.nonzero()
val = X.data

# all products at operation point (equivalent to accumarray in MATLAB)
# Efficiency: vectorized operations with boolean masking, avoids nested loops
Y = np.zeros(r)
for c in range(r):
    mask = coli == c
    if np.any(mask):
        Y[c] = np.prod(val[mask])

# Columns with 2 zero factors will be zero (critical case handling)
# Efficiency: vectorized boolean indexing, O(r) operation
Y[ctwo] = 0

# invert only nonzero elements (efficient: operates only on data array)
X.data = 1.0 / X.data

# compute factor matrix -> Paper
F = S.multiply(Y).multiply(X)

EABC = P * (F.T)  # Compute combined LTI model matrix
print('Linearization: done')  # display when linearization is done
t1 = time.process_time()  # get cputime for linearization

# %% Extract LTI matrices
E = -EABC[:, 0:n]  # Extract "mass" matrix E
A = EABC[:, n:2*n]  # Extract system matrix A
B = EABC[:, 2*n:2*n+m]  # Extract input matrix B

# Local stability: dominant (geralized) eigenvalues
A_dense = A.toarray()
E_dense = E.toarray()

eigvals, eigvecs = sp.linalg.eig(A_dense, E_dense)  # resuelve A v = λ E v
print(eigvals)
toc = time.process_time()
print("Time: ", toc-tic)

