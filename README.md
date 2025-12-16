# Python CPN Linearization for MTI Models

Linearization of MTI models in CPN  representation using sparse matrices.

## Requirements

- **MATLAB**: R2016b or later (with Control System Toolbox)
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

## Running the Scripts

### MATLAB
```matlab
original_gerwald
```

Or from command line:
```bash
matlab -batch "addpath(pwd); original_gerwald"
```

### Python
```bash
python cpn_linearize_basics.py
```

## Configuration (Python)

Edit these flags at the top of `cpn_linearize_basics.py`:

- `test = False`: Use defined matrices (True = use predefined test matrices)
- `test_zeros = False`: Use `dx = ones` (True = use `dx = zeros`, conflictive case)
- `debug = True`: Print intermediate matrices

## What It Does

1. Defines structure matrix `S` and parameter matrix `P`
2. Computes linearization factor matrix `F`
3. Extracts system matrices `E`, `A`, `B`
4. Computes generalized eigenvalues for stability analysis

## Files

- `original_gerwald.m` - Original MATLAB implementation
- `cpn_linearize_basics.py` - Python implementation

## Author

Alexandre Gràcia i Calvo

Original MATLAB: Gerwald Lichtenberg (16.2.2025)  
Python port: TenSyGrid Hackathon 16.12.2025
