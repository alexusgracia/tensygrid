import sys
import os
import numpy as np
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import numpy as np
import sympy as sp
import os
import sys

# Ensure the builder class is accessible
import PolynomialMatrixBuilder_class as PMB

if __name__ == "__main__":
    np.random.seed(1)
    n_vars = 5  # Groups of 5 for x, dx, u, y (20 total variables)
    
    # 1. Generate a 5x5 symmetric matrix H for the state dynamics (A)
    data = np.random.randn(n_vars, n_vars)
    H = (data + data.T) / 2

    # 2. Construct 20 symbolic equations
    eqs = []

    # --- BLOCK 1: Differential Equations (Rows 1-5) ---
    # Target: E*dx = A*x -> Structure: -1.0*dxi + sum(H_ij * xj) = 0
    # The builder multiplies the -1.0 from dx by -1, resulting in +1.0 (Identity in E)
    for i in range(n_vars):
        eq_str = f"-1.0*dx{i+1}"
        for j in range(n_vars):
            coeff = H[i, j]
            eq_str += f" + ({coeff})*x{j+1}"
        eqs.append(eq_str)

    # --- BLOCK 2: Algebraic Equations (Rows 6-20) ---
    # Fills the remaining 15 rows to maintain a 20x20 system
    # Relation y = u (Rows 6-10)
    for i in range(5):
        eqs.append(f"1.0*y{i+1} - 1.0*u{i+1}")
    
    # Relation y = x (Rows 11-20)
    for i in range(10):
        idx = (i % 5) + 1
        eqs.append(f"1.0*y{idx} - 1.0*x{idx}")

    # 3. Initialize the Builder
    builder = PMB.PolynomialMatrixBuilder(eqs, ineqs=[], verbose=False)

    # 4. Define Operating Point (OP)
    v_dict = {s.name: 1.0 for s in builder.all_symbols}
    
    # 5. Execute Linearization
    builder.linearize(v_dict)
    
    # 6. Stability Analysis
    # 5 finite eigenvalues (from H) and 15 infinite eigenvalues (algebraic)
    evals, e_l, e_r, p_m, stable, margin = builder.compute_stability()

    # 7. Final Reporting
    print("\n" + "="*35)
    print("DESCRIPTOR SYSTEM ANALYSIS (20x20)")
    print("="*35)
    
    builder.report(
        eigenvalues=evals,
        is_stable=stable,
        max_real=margin,
        print_matrices=True
    )

    # Validation: Compare finite eigenvalues with theoretical H eigenvalues
    h_evals = np.linalg.eigvals(H)
    print(f"\nTheoretical Eigenvalues of H (should match finite ones):")
    print(np.sort(h_evals))