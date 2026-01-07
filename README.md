# Python CPN Linearization for MTI Models

Linearization of MTI models in CPN  representation using sparse matrices.

## Requirements

- **MATLAB**: R2016b or later (with Control System Toolbox) - for reference implementations
- **Python**: 3.8+ with NumPy and SciPy

## Installation (Python)

1. Create virtual environment:
```bash
python -m venv venv
```

2. Activate it:
   - Windows: `.\venv\Scripts\Activate.ps1`
   - Linux/Mac: `source venv/bin/activate`

3. Install dependencies:
```bash
pip install numpy scipy
```

## Files

### Python Implementations

- **`cpn_linearize_basics.py`** - Basic implementation with core linearization algorithm
  - Handles standard linearization cases
  - Efficient sparse matrix operations
  - Suitable for most use cases
  
- **`cpn_linearize_basics2.py`** - Enhanced implementation with improved edge case handling
  - Explicit computation of critical indices (`izero`, `imone`, `ctwo`)
  - Better alignment with MATLAB reference implementation
  - More robust handling of singular cases
  - Recommended for production use

### MATLAB Reference Implementations

- **`original_files/original_gerwald.m`** - Original MATLAB implementation (basic version)
- **`original_files/second_iteration_gerwald.m`** - Enhanced MATLAB implementation with critical case handling

### Documentation

- **`README_cpn_linearize_basics.md`** - Detailed documentation for basic Python implementation
  - Step-by-step algorithm explanation
  - Function reference with efficiency notes
  - Mathematical background
  
- **`README_cpn_linearize_basics2.md`** - Detailed documentation for enhanced Python implementation
  - Comparison with basic version
  - Advanced function explanations
  - Performance optimizations

## Running the Scripts

### Python

**Basic version:**
```bash
python cpn_linearize_basics.py
```

**Enhanced version:**
```bash
python cpn_linearize_basics2.py
```

### MATLAB

**Original version:**
```matlab
original_gerwald
```

**Enhanced version:**
```matlab
second_iteration_gerwald
```

Or from command line:
```bash
matlab -batch "addpath(pwd); original_gerwald"
matlab -batch "addpath(pwd); second_iteration_gerwald"
```

## Configuration (Python)

### Basic Version (`cpn_linearize_basics.py`)

Edit these flags at the top of the file:

- `check = 0`: 0 = sparse/scalable/fast | 1 = full (only for checks)
- `test = True`: Use predefined test case (True) or random matrices (False)
- `test_zeros = False`: Use `dx = ones` (False) or `dx = zeros` (True, conflictive case)
- `debug = True`: Print intermediate matrices for debugging

### Enhanced Version (`cpn_linearize_basics2.py`)

Same configuration options as basic version, with additional internal optimizations.

## What It Does

The linearization algorithm performs the following steps:

1. **Problem Definition**: Sets up structure matrix `S` and parameter matrix `P` for the iMTI model
2. **Factor Computation**: Computes factor matrix `X = S ⊙ v - |S|` where `v` is the signal vector
3. **Critical Case Handling**: Identifies and corrects special cases (zero factors, singular points)
4. **Product Computation**: Computes products of factors per tensor rank component
5. **Factor Matrix**: Constructs final factor matrix `F` for linearization
6. **Matrix Extraction**: Extracts LTI system matrices:
   - `E`: Mass/descriptor matrix
   - `A`: System matrix
   - `B`: Input matrix
7. **Stability Analysis**: Computes generalized eigenvalues `λ` where `A·v = λ·E·v`

## Key Features

### Sparse Matrix Efficiency
- Operations only on nonzero elements (O(nnz) complexity)
- Scalable to large systems when matrices are sparse
- Memory efficient storage

### Algorithm Variants
- **Basic**: Core algorithm with essential features
- **Enhanced**: Additional robustness for edge cases and critical indices

### Performance
- Sparse operations: O(nnz) where nnz is number of nonzero elements
- Eigenvalue computation: O(n³) for dense matrices (only for small systems)

## Documentation

For detailed explanations, see:

- **[README_cpn_linearize_basics.md](README_cpn_linearize_basics.md)** - Complete guide to basic implementation
  - Step-by-step algorithm walkthrough
  - Function reference (`multiply`, `toarray`, `tocoo`, `tocsr`, etc.)
  - Mathematical background
  - Complexity analysis

- **[README_cpn_linearize_basics2.md](README_cpn_linearize_basics2.md)** - Complete guide to enhanced implementation
  - Differences from basic version
  - Advanced features (`izero`, `imone`, `ctwo`)
  - Performance optimizations
  - Comparison table

## Example Usage

### Basic Linearization

```python
# Set configuration
test = True
debug = False

# Run linearization
python cpn_linearize_basics.py

# Output: Eigenvalues and computation time
```

### Enhanced Linearization with Critical Cases

```python
# Set configuration
test = False
debug = True

# Run enhanced version
python cpn_linearize_basics2.py

# Output: Detailed progress messages, eigenvalues, timing breakdown
```

## Mathematical Background

The algorithm linearizes nonlinear iMTI models of the form:

```
E·dx = f(x, u)
```

where `f(x, u)` is represented in CPN form using structure matrix `S` and parameter matrix `P`. The linearization around operating point `(x₀, u₀)` yields:

```
E·dx = A·x + B·u
```

where `E`, `A`, and `B` are extracted from the factor matrix `F` computed during linearization.

## Authors

- **Python Implementation**: Alexandre Gràcia i Calvo
- **Original MATLAB**: Gerwald Lichtenberg (16.12.2025)
- **Project**: TenSyGrid Hackathon 16.12.2025

## License

[Add license information if applicable]

## References

- CPN (Canonical Polyadic Network) representation theory
- iMTI (implicit Multi-Tensor Index) models
- Sparse matrix linearization techniques
