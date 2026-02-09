# README: cpn_linearize_basics.py

## Overview

This script implements a **basic linearization algorithm for iMTI (implicit Multi-Tensor Index) models in CPN (Canonical Polyadic Network) representation**. The algorithm performs sparse matrix operations to efficiently linearize nonlinear systems around an operating point, extracting Linear Time-Invariant (LTI) system matrices (E, A, B) for stability analysis.

## Table of Contents

1. [Imports and Configuration](#imports-and-configuration)
2. [Problem Definition](#problem-definition)
3. [Linearization Algorithm](#linearization-algorithm)
4. [Matrix Extraction](#matrix-extraction)
5. [Stability Analysis](#stability-analysis)
6. [Key Functions Explained](#key-functions-explained)

---

## Imports and Configuration

### Imports

```python
import numpy as np
import scipy as sp
import time
```

- **`numpy`**: Provides array operations and mathematical functions
- **`scipy`**: Provides sparse matrix operations and linear algebra functions
- **`time`**: Used for performance measurement

### Configuration Variables

```python
check = 0      # 0 = sparse/scalable/fast | 1 = full (only for checks)
test = True    # Use predefined test case when True
test_zeros = False  # Use dx=zeros (conflictive case) when True
debug = True   # Print intermediate matrices when True
```

---

## Problem Definition

### System Dimensions

The script defines an iMTI model with the following dimensions:

- **`n`**: Number of states (state variables)
- **`m`**: Number of inputs (control inputs)
- **`p`**: Number of outputs (currently only 0 supported)
- **`q`**: Number of equations (typically `q = n`)
- **`N`**: Total number of signals = `2*n + m + p`
- **`r`**: Tensor rank (number of columns in structure matrix S)

### Test Case (`test = True`)

When `test = True`, the script uses a predefined simple example:

**System equations:**
- `dx1 = -x1 + u`
- `dx2 = -x2 + x1*x2`

**Structure Matrix (S)**: 
- Shape: `(N, r) = (5, 6)`
- Defines which signals appear in which tensor terms
- Each row corresponds to a signal (dx1, dx2, x1, x2, u)
- Each column corresponds to a tensor rank component

**Parameter Matrix (P)**:
- Shape: `(q, r) = (2, 6)`
- Contains the coefficients for each equation and tensor term
- Row 1: coefficients for equation 1 (dx1)
- Row 2: coefficients for equation 2 (dx2)

**Linearization Point**:
- `dx`: State derivatives (ones or zeros depending on `test_zeros`)
- `x`: States (2*ones)
- `u`: Inputs (ones)
- `y`: Outputs (ones)

### Random Case (`test = False`)

When `test = False`, the script generates random matrices:
- Random structure matrix S with values in {-1, 0, 1}
- Random parameter matrix P (sparse, ~30% density)
- Random linearization point values

### Signal Vector Construction

```python
v = np.vstack([dx, x, u, y])
```

**`np.vstack()`**: Vertically stacks arrays
- Creates a column vector by stacking: state derivatives, states, inputs, outputs
- Result shape: `(N, 1)` where `N = 2*n + m + p`

---

## Linearization Algorithm

### Step 1: Compute Factor Matrix X

```python
X = S.multiply(v) - abs(S)
```

**What this does:**
- Computes `X = S ⊙ v - |S|` where `⊙` is element-wise multiplication
- This represents "all Factors - 1" in the linearization formula

**`S.multiply(v)`**:
- **Purpose**: Element-wise multiplication between sparse matrix S and vector v
- **How it works**: For each nonzero element S[i,j], computes S[i,j] * v[i]
- **Efficiency**: Only operates on nonzero elements (sparse operation)
- **Result**: Sparse matrix of same shape as S

**`abs(S)`**:
- Computes absolute value of each element in S
- Maintains sparsity structure

**Why subtract abs(S)?**
- The linearization formula requires factors of the form `(1 - |S| + S*v)`
- Here we compute `S*v - |S|`, which equals `(1 - |S| + S*v) - 1`
- This "minus 1" will be corrected later

### Step 2: Check for Critical Cases

```python
icrit = X == -1
if icrit.nnz > 0:
    raise RuntimeError("specials to be implemented ?")
```

**`X == -1`**:
- Creates a boolean sparse matrix indicating where X equals -1
- **`.nnz`**: Number of nonzero elements
- If any X[i,j] == -1, this indicates a critical/singular case that needs special handling

### Step 3: Add 1 to All Nonzero Elements

```python
X.data = X.data + 1.0
```

**`X.data`**:
- **Purpose**: Direct access to the data array of a sparse matrix
- **What it contains**: Array of all nonzero values in the sparse matrix
- **Efficiency**: Modifying `.data` directly is O(nnz) instead of O(N*r) for full matrix
- **Result**: All nonzero elements of X are incremented by 1

**Why this works:**
- After `X = S*v - abs(S)`, we need `X = 1 - |S| + S*v`
- Adding 1 to all nonzero elements gives us `(S*v - |S|) + 1 = 1 - |S| + S*v` for nonzero positions

### Step 4: Handle Zero Cases (S != 0 but X == 0)

```python
S_nonzero = S.nonzero()
X_nz_dict = set(zip(*X.nonzero()))
```

**`S.nonzero()`**:
- **Purpose**: Returns indices of all nonzero elements
- **Returns**: Tuple `(row_indices, col_indices)`
- **Example**: If S[0,0]=1 and S[1,2]=-1, returns `([0, 1], [0, 2])`

**`X.nonzero()`**:
- Same as above but for matrix X
- **`zip(*X.nonzero())`**: Transposes the tuple to get pairs `(row, col)`
- **`set(...)`**: Creates a set for O(1) lookup

**The Problem:**
- When `S[i,j] != 0` but `X[i,j] == 0`, we need to set `X[i,j] = 1`
- This happens when `S[i,j]*v[i] == abs(S[i,j])`, making `X[i,j] = 0`
- But X is sparse, so zero elements aren't stored

**The Solution:**
```python
to_add = []
for i, j in zip(*S_nonzero):
    if (i, j) not in X_nz_dict:
        to_add.append((i, j))
```

- Finds all positions where S is nonzero but X is zero (not in X_nz_dict)
- These positions need to be added to X with value 1.0

**Adding New Elements:**
```python
rows, cols = zip(*to_add)
data = np.ones(len(to_add))
X_add = sp.sparse.coo_matrix((data, (rows, cols)), shape=shape)
X = X + X_add
```

**`sp.sparse.coo_matrix()`**:
- **Purpose**: Creates a sparse matrix in COO (Coordinate) format
- **Format**: Stores (row, col, value) triplets
- **Arguments**: `(data, (rows, cols), shape)`
  - `data`: Array of values
  - `(rows, cols)`: Tuple of row and column index arrays
  - `shape`: Matrix dimensions
- **Why COO**: Easy to construct from coordinate lists

**`X = X + X_add`**:
- Adds the new elements to X
- Automatically handles duplicate positions

**`X.tocsr()`**:
- **Purpose**: Converts sparse matrix to CSR (Compressed Sparse Row) format
- **CSR Format**: Efficient for row operations and matrix-vector products
- **Structure**: 
  - `data`: Array of nonzero values
  - `indices`: Column indices
  - `indptr`: Row pointer array (indicates where each row starts)
- **Why convert**: CSR is more efficient for subsequent operations

### Step 5: Compute Products Y

```python
rowi, coli = X.nonzero()
val = X.data
Y = np.zeros(r)
for c in range(r):
    mask = coli == c
    if np.any(mask):
        Y[c] = np.prod(val[mask])
```

**What this does:**
- Computes the product of all nonzero elements in each column of X
- This is equivalent to MATLAB's `accumarray(coli, val, [r 1], @prod)`

**Step-by-step:**
1. **`X.nonzero()`**: Get all (row, col) indices of nonzero elements
2. **`X.data`**: Get all nonzero values
3. **For each column c**:
   - **`mask = coli == c`**: Boolean array indicating which elements are in column c
   - **`np.any(mask)`**: Check if column c has any nonzero elements
   - **`np.prod(val[mask])`**: Product of all values in column c

**Why this matters:**
- Y[c] represents the product of all factors for tensor rank component c
- This is a key step in the CPN linearization formula

### Step 6: Invert X

```python
X.data = 1.0 / X.data
```

**Purpose**: Invert all nonzero elements of X
- **Efficiency**: Direct operation on data array, O(nnz) complexity
- **Result**: X[i,j] becomes 1/X[i,j] for all nonzero positions

### Step 7: Compute Factor Matrix F

```python
F = S.multiply(Y).multiply(X)
```

**What this does:**
- Computes the final factor matrix F
- **`S.multiply(Y)`**: Element-wise multiply S by vector Y (broadcasting)
  - For each S[i,j], computes S[i,j] * Y[j]
- **`.multiply(X)`**: Element-wise multiply result by X
  - Final result: F[i,j] = S[i,j] * Y[j] * X[i,j]

**Mathematical meaning:**
- F[i,j] represents the linearization coefficient for signal i in tensor term j
- This is the derivative of the nonlinear function with respect to signal i

---

## Matrix Extraction

### Step 8: Compute Combined LTI Matrix

```python
EABC = P * (F.T)
```

**`F.T`**:
- **Purpose**: Transpose of matrix F
- **Why transpose**: P has shape (q, r), F has shape (N, r)
- We need to multiply P with F^T to get shape (q, N)

**Matrix multiplication**:
- **`P * F.T`**: Standard matrix multiplication
- **Result shape**: (q, N) = (n, 2*n + m + p)
- Contains all LTI matrices concatenated: [E, A, B]

### Step 9: Extract Individual Matrices

```python
E = -EABC[:, 0:n]        # Mass matrix (descriptor system)
A = EABC[:, n:2*n]      # System matrix
B = EABC[:, 2*n:2*n+m]  # Input matrix
```

**Matrix slicing**:
- **`EABC[:, 0:n]`**: First n columns (correspond to state derivatives dx)
- **`EABC[:, n:2*n]`**: Next n columns (correspond to states x)
- **`EABC[:, 2*n:2*n+m]`**: Next m columns (correspond to inputs u)

**Why negative for E?**
- The descriptor system form is: `E*dx = A*x + B*u`
- Rearranging: `-E*dx + A*x + B*u = 0`
- The sign convention requires E to be negated

---

## Stability Analysis

### Step 10: Compute Eigenvalues

```python
A_dense = A.toarray()
E_dense = E.toarray()
eigvals, eigvecs = sp.linalg.eig(A_dense, E_dense)
```

**`A.toarray()` / `E.toarray()`**:
- **Purpose**: Converts sparse matrix to dense (full) NumPy array
- **Why needed**: `sp.linalg.eig()` requires dense matrices for generalized eigenvalue problems
- **Trade-off**: Memory vs. computation speed
- **When to use**: Small matrices (n < 1000) or when dense operations are required

**`sp.linalg.eig(A, E)`**:
- **Purpose**: Solves generalized eigenvalue problem `A*v = λ*E*v`
- **Returns**:
  - `eigvals`: Array of eigenvalues (complex numbers)
  - `eigvecs`: Matrix of eigenvectors (columns)
- **Stability criterion**: System is stable if all eigenvalues have negative real parts

**Generalized Eigenvalue Problem:**
- Standard form: `A*v = λ*v`
- Generalized form: `A*v = λ*E*v`
- When E = I (identity), reduces to standard form
- E is called the "mass matrix" or "descriptor matrix"

---

## Key Functions Explained

### Sparse Matrix Operations

#### `multiply(other)`
- **Purpose**: Element-wise multiplication
- **Syntax**: `A.multiply(B)` or `A * B` (if B is scalar/vector)
- **Efficiency**: O(nnz) - only operates on nonzero elements
- **Example**: If A[0,0]=2 and B[0,0]=3, result[0,0]=6

#### `toarray()`
- **Purpose**: Convert sparse matrix to dense NumPy array
- **Memory**: Creates full N×M array (even for zeros)
- **Use case**: When dense operations are required
- **Trade-off**: Memory usage vs. operation compatibility

#### `tocoo()`
- **Purpose**: Convert to COO (Coordinate) format
- **Format**: Stores (row, col, value) triplets
- **Advantages**: Easy to construct, modify, and convert from
- **Disadvantages**: Slower for arithmetic operations
- **Use case**: Building matrices from coordinate lists

#### `tocsr()`
- **Purpose**: Convert to CSR (Compressed Sparse Row) format
- **Format**: Row-compressed storage
- **Advantages**: Fast row operations, matrix-vector products
- **Use case**: Final format for computation

#### `nonzero()`
- **Purpose**: Get indices of nonzero elements
- **Returns**: Tuple `(row_indices, col_indices)`
- **Efficiency**: O(nnz) - only iterates over nonzero elements
- **Use case**: Finding positions to modify or analyze

#### `.data` attribute
- **Purpose**: Direct access to nonzero values array
- **Efficiency**: O(1) access, O(nnz) modification
- **Use case**: Bulk operations on nonzero values
- **Warning**: Modifying `.data` directly can break matrix invariants if not careful

### NumPy Functions

#### `np.vstack(arrays)`
- **Purpose**: Vertically stack arrays
- **Example**: `vstack([a, b])` creates `[a; b]` (MATLAB notation)
- **Use case**: Concatenating vectors into single column vector

#### `np.bincount(x, minlength)`
- **Purpose**: Count occurrences of each value in array
- **Returns**: Array where result[i] = count of value i in x
- **`minlength`**: Ensures output has at least this length
- **Use case**: Counting how many times each column index appears

#### `np.prod(array)`
- **Purpose**: Product of all elements in array
- **Example**: `prod([2, 3, 4]) = 24`
- **Use case**: Computing products for tensor terms

#### `np.ones(shape)`
- **Purpose**: Create array filled with ones
- **Example**: `ones((3, 1))` creates 3×1 array of ones
- **Use case**: Initializing vectors with unit values

#### `np.zeros(shape)`
- **Purpose**: Create array filled with zeros
- **Example**: `zeros(r)` creates array of length r
- **Use case**: Initializing result arrays

#### `np.any(array)`
- **Purpose**: Check if any element is True/nonzero
- **Returns**: Boolean
- **Use case**: Checking if mask has any matches

---

## Algorithm Complexity

- **Space**: O(nnz) where nnz is number of nonzero elements
- **Time**: 
  - Sparse operations: O(nnz)
  - Dense conversions: O(N×r) when needed
  - Eigenvalue computation: O(n³) for dense matrices

## Notes

- This is the **basic version** - it handles the core linearization algorithm
- The algorithm is **scalable** for sparse systems (works with large n when matrices are sparse)
- **Critical cases** (X == -1) are detected but not fully handled (raises error)
- The **zero case handling** (S != 0 but X == 0) is implemented efficiently using set operations

