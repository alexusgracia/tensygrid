# README: cpn_linearize_basics2.py

## Overview

This script implements an **enhanced linearization algorithm for iMTI (implicit Multi-Tensor Index) models in CPN (Canonical Polyadic Network) representation**. It extends the basic version (`cpn_linearize_basics.py`) with improved handling of critical cases, including explicit computation of `izero`, `imone`, and `ctwo` indices for more robust linearization. This version closely follows the MATLAB implementation in `second_iteration_gerwald.m`.

## Table of Contents

1. [Key Differences from Basic Version](#key-differences-from-basic-version)
2. [Imports and Configuration](#imports-and-configuration)
3. [Problem Definition](#problem-definition)
4. [Enhanced Linearization Algorithm](#enhanced-linearization-algorithm)
5. [Matrix Extraction](#matrix-extraction)
6. [Stability Analysis](#stability-analysis)
7. [Advanced Functions Explained](#advanced-functions-explained)

---

## Key Differences from Basic Version

### Main Enhancements

1. **Explicit `izero` computation**: Identifies positions where `S.*v == abs(S)` and `S != 0`
2. **Explicit `imone` computation**: Identifies positions where `X == -1` (before adding 1)
3. **`ctwo` computation**: Identifies columns with more than one `X == -1`
4. **Combined correction**: Handles both `izero` and `imone` cases in a single operation
5. **Critical column handling**: Sets `Y[ctwo] = 0` for columns with multiple zero factors
6. **Performance comments**: Detailed efficiency explanations for each optimization

### Algorithmic Improvements

- More robust handling of edge cases
- Better alignment with MATLAB reference implementation
- Explicit tracking of critical indices for debugging and analysis

---

## Imports and Configuration

### Imports

```python
import numpy as np
import scipy as sp
import time
```

Same as basic version - see [README_cpn_linearize_basics.md](README_cpn_linearize_basics.md) for details.

### Configuration Variables

```python
check = 0      # 0 = sparse/scalable/fast | 1 = full (only for checks)
test = False  # Use predefined test case when True
test_zeros = False  # Use dx=zeros (conflictive case) when True
debug = False # Print intermediate matrices when True
```

---

## Problem Definition

Identical to basic version - see [README_cpn_linearize_basics.md](README_cpn_linearize_basics.md) for complete details.

---

## Enhanced Linearization Algorithm

### Step 1: Compute Factor Matrix X

```python
X = S.multiply(v) - abs(S)  # sparse matrix of all Factors - 1
```

Same as basic version - computes `X = S ⊙ v - |S|`.

### Step 2: Compute Critical Indices (`izero`, `imone`, `ctwo`)

This is the **key enhancement** over the basic version. The algorithm explicitly computes three types of critical indices:

#### 2.1: Compute `izero` (Indices where S.*v == abs(S) and S != 0)

```python
Sv = S.multiply(v)
absS = abs(S)
Sv_dense = Sv.toarray()
absS_dense = absS.toarray()
S_dense = S.toarray()
izero_mask = (np.abs(Sv_dense - absS_dense) < 1e-10) & (S_dense != 0)
izero = sp.sparse.csr_matrix(izero_mask.astype(float))
```

**What `izero` represents:**
- Positions where `S[i,j]*v[i] == abs(S[i,j])`
- This means `X[i,j] = S[i,j]*v[i] - abs(S[i,j]) = 0`
- But `S[i,j] != 0`, so this is a position where S is nonzero but X became zero
- These positions need to be set to 1 in X

**Why use dense comparison:**
- **Accuracy**: Sparse comparison (`Sv == absS`) can miss zero values
- **Completeness**: Dense comparison checks all positions, not just stored nonzero elements
- **Trade-off**: Memory vs. accuracy - acceptable for small to medium matrices

**`np.abs(Sv_dense - absS_dense) < 1e-10`:**
- **Purpose**: Numerical comparison with tolerance
- **Why tolerance**: Floating-point arithmetic may not give exact equality
- **1e-10**: Typical tolerance for numerical comparisons

**`& (S_dense != 0)`:**
- **Purpose**: Ensure we only consider positions where S is nonzero
- **Why**: We only care about positions where S has structure

**`sp.sparse.csr_matrix(izero_mask.astype(float))`:**
- Converts boolean mask to sparse matrix
- Only stores True positions (as 1.0)
- CSR format for efficient operations

#### 2.2: Compute `imone` (Indices where X == -1)

```python
imone = sp.sparse.csr_matrix((X == -1).astype(float))
```

**What `imone` represents:**
- Positions where `X[i,j] == -1` **before** adding 1
- This happens when `S[i,j]*v[i] - abs(S[i,j]) == -1`
- After adding 1, these become 0, which is problematic
- These positions also need to be set to 1

**Why sparse boolean operation:**
- **Efficiency**: Only stores positions where condition is True
- **Memory**: O(nnz_critical) instead of O(N×r)
- **Speed**: Fast comparison on sparse structure

#### 2.3: Compute `ctwo` (Columns with more than one X == -1)

```python
imone_rows, imone_cols = imone.nonzero()
if len(imone_cols) > 0:
    col_counts = np.bincount(imone_cols, minlength=r)
    ctwo = col_counts > 1
else:
    ctwo = np.zeros(r, dtype=bool)
```

**What `ctwo` represents:**
- Boolean array of length r
- `ctwo[j] = True` if column j has **more than one** `X == -1`
- These columns will have zero product in Y (multiple zero factors)

**`imone.nonzero()`:**
- Gets all (row, col) positions where `X == -1`
- Returns tuple `(row_indices, col_indices)`

**`np.bincount(imone_cols, minlength=r)`:**
- **Purpose**: Count how many times each column index appears
- **How it works**: 
  - Creates array of size `max(imone_cols) + 1` (or `minlength`)
  - `result[i]` = number of times value `i` appears in `imone_cols`
- **`minlength=r`**: Ensures output has at least r elements (pads with zeros)
- **Example**: If `imone_cols = [0, 0, 2, 2, 2]`, then `bincount` returns `[2, 0, 3, ...]`

**`col_counts > 1`:**
- Creates boolean array where True indicates columns with count > 1
- These are the problematic columns

**Why this matters:**
- If a column has multiple `X == -1`, after adding 1 they become 0
- Product of factors including zeros is zero
- We need to explicitly set `Y[ctwo] = 0` later

### Step 3: Check for Critical Cases

```python
icrit = X == -1  # get numerical critical indices
if icrit.nnz > 0:
    raise RuntimeError("specials to be implemented ?")
```

Same as basic version - detects if any critical cases exist.

### Step 4: Add 1 to All Nonzero Elements

```python
X.data = X.data + 1.0
```

Same as basic version - increments all nonzero elements.

### Step 5: Combined Correction (izero and imone)

This is a **key improvement** - handles both correction cases in one operation:

```python
izero_coo = izero.tocoo()
imone_coo = imone.tocoo()
X_coo = X.tocoo()
```

**Why convert to COO:**
- **COO format**: Easy to access and modify individual elements
- **Row/Col access**: Direct access to `.row`, `.col`, `.data` arrays
- **Modification**: Can append new elements easily

**Get all correction positions:**
```python
izero_positions = set(zip(izero_coo.row, izero_coo.col))
imone_positions = set(zip(imone_coo.row, imone_coo.col))
all_correction_positions = izero_positions | imone_positions
```

**`zip(izero_coo.row, izero_coo.col)`:**
- Pairs row and column indices: `(row0, col0), (row1, col1), ...`
- Creates coordinate pairs for each nonzero element

**`set(...)`:**
- **Purpose**: Creates set for O(1) membership testing
- **Efficiency**: Set union `|` is O(n) where n is number of elements
- **Why sets**: Avoid duplicates when izero and imone overlap

**`izero_positions | imone_positions`:**
- **Purpose**: Set union - all positions that need correction
- **Efficiency**: O(n) where n is total unique positions

**Apply corrections:**
```python
if all_correction_positions:
    X_positions = {}
    for idx, (row_idx, col_idx) in enumerate(zip(X_coo.row, X_coo.col)):
        X_positions[(row_idx, col_idx)] = idx
```

**Create position mapping:**
- **Purpose**: Map (row, col) pairs to indices in X_coo.data array
- **Efficiency**: Dictionary lookup O(1) vs. linear search O(n)
- **Why needed**: To quickly find which data element to modify

**Update existing positions:**
```python
for (i, j) in all_correction_positions:
    if (i, j) in X_positions:
        X_coo.data[X_positions[(i, j)]] = 1.0
    else:
        new_rows.append(i)
        new_cols.append(j)
        new_data.append(1.0)
```

**Two cases:**
1. **Position exists in X**: Modify existing data value to 1.0
2. **Position doesn't exist**: Need to add new element (collect for batch addition)

**Add new positions:**
```python
if new_rows:
    X_coo.row = np.concatenate([X_coo.row, new_rows])
    X_coo.col = np.concatenate([X_coo.col, new_cols])
    X_coo.data = np.concatenate([X_coo.data, new_data])
    X = X_coo.tocsr()
```

**`np.concatenate()`:**
- **Purpose**: Concatenate arrays along specified axis
- **Efficiency**: Single operation instead of multiple appends
- **Why batch**: More efficient than adding elements one by one

**Convert back to CSR:**
- CSR format is more efficient for subsequent operations
- COO is good for construction, CSR is good for computation

### Step 6: Compute Products Y

```python
rowi, coli = X.nonzero()
val = X.data
Y = np.zeros(r)
for c in range(r):
    mask = coli == c
    if np.any(mask):
        Y[c] = np.prod(val[mask])
```

Same as basic version - computes product of factors per column.

### Step 7: Handle Critical Columns (`ctwo`)

```python
Y[ctwo] = 0
```

**This is NEW** - explicitly handles columns with multiple zero factors:

**What this does:**
- Sets `Y[j] = 0` for all columns j where `ctwo[j] = True`
- These columns had multiple `X == -1`, which became zeros after adding 1
- Product of factors including zeros is zero

**Why explicit:**
- Even if the product computation gives zero, explicitly setting it ensures correctness
- Handles edge cases where numerical precision might not give exact zero
- Makes the algorithm more robust

**`Y[ctwo]`:**
- **Boolean indexing**: Selects elements where `ctwo` is True
- **Efficiency**: O(r) operation, very fast
- **Vectorized**: NumPy handles this efficiently

### Step 8: Invert X

```python
X.data = 1.0 / X.data
```

Same as basic version - inverts all nonzero elements.

### Step 9: Compute Factor Matrix F

```python
F = S.multiply(Y).multiply(X)
```

Same as basic version - computes final factor matrix.

---

## Matrix Extraction

### Step 10: Compute Combined LTI Matrix

```python
EABC = P * (F.T)  # Compute combined LTI model matrix
print('Linearization: done')  # display when linearization is done
t1 = time.process_time()  # get cputime for linearization
```

**New additions:**
- **Progress message**: Indicates when linearization completes
- **Timing**: Records time for performance analysis
- **`t1`**: Used later to report timing breakdown

### Step 11: Extract Individual Matrices

```python
# %% Extract LTI matrices
E = -EABC[:, 0:n]  # Extract "mass" matrix E
A = EABC[:, n:2*n]  # Extract system matrix A
B = EABC[:, 2*n:2*n+m]  # Extract input matrix B
```

Same as basic version - extracts E, A, B matrices.

---

## Stability Analysis

### Step 12: Compute Eigenvalues

```python
A_dense = A.toarray()
E_dense = E.toarray()
eigvals, eigvecs = sp.linalg.eig(A_dense, E_dense)  # resuelve A v = λ E v
print(eigvals)
toc = time.process_time()
print("Time: ", toc-tic)
```

Same as basic version - computes generalized eigenvalues for stability analysis.

---

## Advanced Functions Explained

### Additional Functions Beyond Basic Version

#### `np.bincount(x, minlength)`

**Purpose**: Count occurrences of each value in an array

**How it works:**
```python
x = [0, 0, 2, 2, 2, 5]
counts = np.bincount(x, minlength=6)
# Result: [2, 0, 3, 0, 0, 1]
# Meaning: value 0 appears 2 times, value 2 appears 3 times, etc.
```

**Parameters:**
- `x`: Input array (must be non-negative integers)
- `minlength`: Minimum size of output array (pads with zeros)

**Use case in this code:**
- Count how many times each column index appears in `imone_cols`
- Determines which columns have multiple `X == -1`

**Efficiency:**
- O(n) where n is length of input array
- Very fast for integer arrays

#### `np.concatenate(arrays)`

**Purpose**: Join arrays along specified axis

**Example:**
```python
a = [1, 2, 3]
b = [4, 5]
result = np.concatenate([a, b])  # [1, 2, 3, 4, 5]
```

**Use case in this code:**
- Batch addition of new elements to COO matrix
- More efficient than appending one by one

**Efficiency:**
- O(n) where n is total elements
- Single allocation vs. multiple reallocations

#### Set Operations

**`set()` creation:**
- **Purpose**: Create unordered collection of unique elements
- **Properties**: O(1) average-case membership testing
- **Use case**: Fast lookup of coordinate pairs

**Set union `|`:**
- **Purpose**: Combine two sets (all unique elements from both)
- **Efficiency**: O(n + m) where n, m are set sizes
- **Use case**: Combine izero and imone positions

**Why sets for coordinates:**
- Coordinate pairs `(i, j)` are hashable (can be set elements)
- Fast to check if position exists
- Automatic duplicate removal

#### Dictionary for Position Mapping

**`X_positions = {}`:**
- **Purpose**: Map (row, col) pairs to data array indices
- **Key**: `(row, col)` tuple
- **Value**: Index in `X_coo.data` array

**Why dictionary:**
- **O(1) lookup**: Find data index for given position instantly
- **Alternative**: Linear search would be O(nnz)

**Example:**
```python
X_positions = {(0, 1): 0, (1, 2): 1, (2, 0): 2}
# Position (1, 2) corresponds to X_coo.data[1]
```

---

## Algorithm Complexity Analysis

### Space Complexity
- **Sparse matrices**: O(nnz) where nnz is number of nonzero elements
- **Dense conversions**: O(N×r) when needed (for izero computation)
- **Auxiliary structures**: O(nnz_critical) for izero, imone, ctwo

### Time Complexity
- **Sparse operations**: O(nnz)
- **Dense conversions**: O(N×r) - only for izero computation
- **Set operations**: O(nnz_critical) for position unions
- **Dictionary operations**: O(nnz) for mapping creation, O(1) per lookup
- **Product computation**: O(nnz) - one pass through columns
- **Eigenvalue computation**: O(n³) for dense matrices

### Performance Optimizations

1. **Sparse operations**: Only operate on nonzero elements
2. **Set-based lookup**: O(1) position checking
3. **Dictionary mapping**: O(1) data index lookup
4. **Batch concatenation**: Single operation instead of multiple appends
5. **Vectorized operations**: Boolean indexing, bincount, etc.

---

## Comparison with Basic Version

| Feature | Basic Version | Enhanced Version (basics2) |
|---------|--------------|---------------------------|
| `izero` computation | Implicit (via set difference) | Explicit (dense comparison) |
| `imone` computation | Implicit (via icrit check) | Explicit (sparse boolean) |
| `ctwo` computation | Not computed | Explicit (bincount) |
| Critical column handling | Not explicit | `Y[ctwo] = 0` |
| Correction method | Set difference approach | Combined izero\|imone |
| Performance comments | Minimal | Extensive |
| MATLAB alignment | Partial | Close alignment |

---

## Notes

- This version is **more robust** for edge cases
- **Better alignment** with MATLAB reference implementation
- **Slightly more memory** usage due to dense conversions for izero
- **More explicit** about handling critical cases
- **Better documented** with efficiency explanations
- Suitable for **production use** where robustness is important

---

## References

- MATLAB reference: `second_iteration_gerwald.m`
- Basic version: `cpn_linearize_basics.py`
- CPN representation theory and iMTI models

