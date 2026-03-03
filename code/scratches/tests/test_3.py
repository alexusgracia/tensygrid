import sys
import os

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import PolynomialMatrixBuilder_class as PMB

if __name__ == "__main__":
    import os
    save_output = False   # ← set False to skip txt export

    eqs=['-xp1 + (x3+ 0.2732*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1))+2*pi*50)*(-y2)',
     '-xp2 + (x3+ 0.2732*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1))+2*pi*50)*y1',
     '-xp3 + (0.022508948*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1)))',   
     'y1 - x1',
     'y2- x2']



    builder = PMB.PolynomialMatrixBuilder(eqs, ineqs=[], verbose=False)

    xpOP = [0, 0, 0]
    xOP  = [-0.5737, -0.8193, -0.1597]
    uOP  = [-100.6109, -217.5716, 318.1824]
    yOP  = [-0.5737, -0.8193]

    v_dict = {
        **{f'dx{i}': xpOP[i - 1] for i in range(1, len(xpOP) + 1)},
        **{f'x{i}':  xOP[i - 1]  for i in range(1, len(xOP) + 1)},
        **{f'u{i}':  uOP[i - 1]  for i in range(1, len(uOP) + 1)},
        **{f'y{i}':  yOP[i - 1]  for i in range(1, len(yOP) + 1)},
        **{f'xp{i}': xpOP[i - 1] for i in range(1, len(xpOP) + 1)},
    }


    builder.linearize(v_dict)

    evals, e_left, e_right, p_matrix, stable, margin = builder.compute_stability()

    out_path = os.path.join(os.path.dirname(__file__), "matrices_output.txt") if save_output else None
    builder.report(
        eigenvalues    = evals,
        is_stable      = stable,
        max_real       = margin,
        print_matrices = True,
        save_path      = out_path,
    )
