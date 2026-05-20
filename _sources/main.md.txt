# TenSyGrid

Python toolkit for symbolic linearization of nonlinear dynamical systems using the **iMTI / CPN (Canonical Polynomial Network)** representation.

Given a set of polynomial differential equations, TenSyGrid automatically builds the structural matrices **S** (or S_H / S_W) and **Phi** (or P / Phi_H / Phi_W) used in the CPN formalism, and computes the analytic linearization matrices **E**, **A**, **B**.

---

## Background

The CPN method represents nonlinear polynomial systems in the form:

$$
S \cdot \dot{z} = P \cdot z
$$

where $z$ is a vector collecting all monomials that appear in the equations (states, derivatives, inputs, and higher-order product terms). The matrices $S$ and $P$ (or $\Phi$) encode the structure of the system and allow systematic analytic linearization around an operating point.

For example, the system:

$$
\dot{x}_1 = -x_1 + u, \quad \dot{x}_2 = -x_2 + x_1 x_2
$$

produces the augmented monomial vector $[\dot{x}_1,\, \dot{x}_2,\, x_1,\, x_2,\, u,\, x_1 x_2]$ and the corresponding $S$ / $P$ matrices automatically.

---

## Repository Structure

```
tensygrid/
├── code/
│   ├── main.py                          # Entry point and usage examples
│   ├── MatrixBuilder.py                 # High-level MatrixBuilder class
│   └── scratches/
│       ├── PolynomialMatrixBuilder_class.py   # Core implementation (sparse matrices)
│       ├── tests.py                     # Unit tests for the builder
│       └── notebook_example.ipynb       # Interactive example notebook
├── manual/                              # Sphinx documentation source
├── MTI-toolbox-2-1/                     # Reference MATLAB MTI toolbox
├── apunts.md                            # Development notes (Catalan)
└── old_old/                             # Legacy Python implementations
```

---

## Installation

### Requirements

- Python 3.10+
- [SymPy](https://www.sympy.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)

### Setup

```bash
# Clone the repository
git clone <repo-url>
cd tensygrid

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate        # Linux / macOS
# .\venv\Scripts\Activate.ps1   # Windows PowerShell

# Install dependencies
pip install sympy numpy scipy
```

---

## Quick Start

Equations are passed as strings with the convention that **`dx<n>`** denotes $\dot{x}_n$:

```python
from code.scratches.PolynomialMatrixBuilder_class import PolynomialMatrixBuilder

eqs = [
    "-dx1 - x1 + u",
    "-dx2 - x2 + x1*x2",
]
ineqs = []   # inequality constraints (optional)

builder = PolynomialMatrixBuilder(eqs, ineqs, verbose=True)

print("S_H =")
print(builder.S_H)

print("Phi_H =")
print(builder.Phi_H)
```

### Linearization

After building the matrices, set an operating point and compute the linearized **E**, **A**, **B** matrices:

```python
op = {'dx1': 0.0, 'x1': 1.0, 'u': 1.0, 'dx2': 0.0, 'x2': 1.0}
result = builder.linearize(op)
print(result)   # {'E': ..., 'A': ..., 'B': ...}
```

---

## Features

- **Symbolic parsing** of polynomial equations via SymPy (string or LaTeX input).
- **Automatic monomial extraction** — identifies all product terms and adds them to the augmented state vector.
- **Sparse matrix construction** (SciPy CSC / CSR) for efficiency with large systems.
- **Equality and inequality constraints** (S_H / Phi_H for equalities, S_W / Phi_W for inequalities).
- **Analytic linearization** computing the descriptor system matrices E, A, B at a given operating point.

---

## Reference MATLAB Toolbox

The `MTI-toolbox-2-1/` folder contains the original MATLAB MTI toolbox used as the reference implementation. The Python code in this repository is a reimplementation of its core linearization routines.


