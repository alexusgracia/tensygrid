import numpy as np
import scipy as sp
from scipy import sparse
from scipy.linalg import eig
import sys
import time
import control as ct


# %% Start file for TenSyGrid Hackathon 16.12.2025
# Linearization for iMTI models in CPN representation (basic, functional)

# Configuration
check = 0  # 0 = sparse/scalable/fast | 1 = full (only for checks)

# %% Define problem

# Dimensions
n = 2               # number of states
m = 1               # number of inputs
p = 0               # number of outputs

q = n               # number of equations
N = 2 * n + m + p   # total number of signals

# iMTI model in CPN representation

S = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 1],
    [0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 1, 0],
]) # Structure matrix

r = S.shape[1] # rank (number of non-zero elements in the structure matrix)
S= sparse.csr_matrix(S)

P = np.array([
    [-1, 0, -1, 0, 1, 0],
    [0, -1, 0, -1, 0, 1],
]) # Parameter matrix

P = sparse.csr_matrix(P)

# linearization point
# dx = np.ones((n, 1))                 # state derivative 
dx = np.zeros((n, 1))                 # state derivative 
x = 2*np.ones((n, 1))                  # state
u = np.ones((m, 1))                  # input
y = np.ones((p, 1))                  # output
v = np.vstack([dx, x, u, y])         # signal vector


tic = time.process_time()   # start the cputime clock
# %% Linearization (Scalable = sparse & low rank)
print(S.toarray())
X = S.multiply(v) - abs(S)
print(X.toarray())
icrit = X == -1 # get numerical critical indices
if icrit.nnz > 0:
    raise RuntimeError("specials to be implemented ?")

X.data = X.data + 1.0 # add one only to the nonzero elements


# puedes hacer un codigo para que en todas las posiciones de S que sean diferentes de 0 y en X sean 0 haya un 1, pero no queremos recorrer todas las posicones, porque es una sparse matrix. hazlo de la forma más eficiente posible
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
    X_add = sparse.coo_matrix((data, (rows, cols)), shape=shape)
    X = X + X_add
    X = X.tocsr()  # ensure CSR format for downstream code

print(X.toarray())

# get the indices and values of X
rowi, coli = X.nonzero()
val = X.data
# print(rowi);
# print(coli);
# print(val);

Y = np.zeros(r) # initialize product vector
for c in range(r): # product over nonzero elements
    mask = coli == c # get the indices of the non-zero elements in the column c
    if np.any(mask): # if there are non-zero elements in the column c
        Y[c] = np.prod(val[mask]) # compute the product of the non-zero elements in the column c
# print(Y);

X.data = 1.0 / X.data # Invert only nonzero elements of X
# print(X.toarray())
F = S.multiply(Y).multiply(X)
# print(F.toarray())
# print(F.T.toarray())

EABC = P*(F.T) 
print(EABC.toarray())# Compute combined matrix
E = -EABC[:, 0:n] # state matrix
# print(E.toarray())
A = EABC[:, n:2*n] # system matrix
# print(A.toarray())
B = EABC[:, 2*n:2*n+m] # input matrix
# print(B.toarray())

# Local stability: dominant (geralized) eigenvalues
A_dense = A.toarray()
E_dense = E.toarray()

eigvals, eigvecs = eig(A_dense, E_dense)  # resuelve A v = λ E v
print(eigvals)
