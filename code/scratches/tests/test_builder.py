from code.scratches.PolynomialMatrixBuilder_class import PolynomialMatrixBuilder

import pytest
import numpy as np
from scipy import sparse

# =====================================================================
# HELPER FUNCTION TO REBUILD MATRICES FROM CONSOLE OUTPUT
# =====================================================================
def build_matrix(shape, coords_and_values):
    """Rebuilds a dense matrix from (row, column, value) tuples."""
    if not coords_and_values:
        return None
    rows = [cv[0] for cv in coords_and_values]
    cols = [cv[1] for cv in coords_and_values]
    data = [cv[2] for cv in coords_and_values]
    return sparse.coo_matrix((data, (rows, cols)), shape=shape).toarray()

# =====================================================================
# EXAMPLE 1 DATA (Nonlinear system with fractions)
# =====================================================================
eqs_1 = [
    "3*dx1*y1 + 6*(1/2+1/2*z1)*(2/3-1/3*u1)*x1",
    "y1 - x1"
]
v_dict_1 = {'dx1': 1.0, 'x1': 2.0, 'u1': 3.0, 'y1': 4.0, 'z1': 5.0}

expected_E_1 = build_matrix((2, 3), [
    (0, 0, -12.0)
])

expected_A_1 = build_matrix((2, 3), [
    (0, 0, -5.999999999999998), (1, 0, -1.0),
    (0, 1, 3.0), (1, 1, 1.0),
    (0, 2, -1.9999999999999996)
])
expected_max_real_1 = None # The system does not calculate stability (None)

# =====================================================================
# EXAMPLE 2 DATA (Large 20-state system)
# =====================================================================
xOP_2 = [-0.1668, -0.4590,  0.5462, -1.3023, -0.4115, -0.5401,  1.3665,  0.7068,  0.6076,  0.9198, 
          0.8994,  0.0327, -1.4159, -1.6809,  0.0378,  1.6635,  0.8310,  1.0477,  1.7938, -0.4851]
uOP_2 = [-0.6336,  1.1019, -0.7170,  2.2000,  0.8906, -0.6354]

# TODO: Insert your actual equations for Example 2 here
eqs_2 = ["dx1 - x1"] 

v_dict_2 = {
    **{f'dx{i}': 0.0        for i in range(1, 21)},
    **{f'x{i}':  xOP_2[i - 1]  for i in range(1, 21)},
    **{f'u{i}':  uOP_2[i - 1]  for i in range(1, 7)},
    **{f'xp{i}': 0.0        for i in range(1, 21)}
}

expected_E_2 = build_matrix((20, 20), [
    (0, 0, 1.0), (9, 1, 1.0), (10, 2, 1.0), (11, 3, 1.0), (12, 4, 1.0),
    (13, 5, 1.0), (14, 6, 1.0), (15, 7, 1.0), (16, 8, 1.0), (17, 9, 1.0),
    (18, 10, 1.0), (1, 11, 1.0), (19, 12, 1.0), (2, 13, 1.0), (3, 14, 1.0),
    (4, 15, 1.0), (5, 16, 1.0), (6, 17, 1.0), (7, 18, 1.0), (8, 19, 1.0)
])

expected_A_2 = build_matrix((20, 20), [
    (9, 0, 0.0276), (15, 0, 0.1236), (0, 1, 0.344), (1, 1, 0.6549), 
    (7, 1, 0.0515), (6, 2, 0.2301), (2, 3, 0.0651), (6, 3, 0.2002), 
    (8, 3, 0.2754), (10, 3, 0.2733), (18, 3, 0.4517), (1, 4, 0.7845), 
    (12, 4, 0.0839), (1, 5, 0.004), (4, 6, 0.38984), (14, 6, 0.4065), 
    (2, 7, 0.3473), (5, 7, 0.0856), (10, 7, 0.0994), (3, 8, 0.2887), 
    (3, 9, 0.0327), (4, 9, 0.1381), (5, 9, 0.2286), (15, 9, 0.1613), 
    (1, 10, 0.2282), (9, 10, 0.2263), (1, 11, 0.0038), (17, 11, 0.0578), 
    (10, 12, 0.1905), (12, 12, 0.2854), (10, 13, 0.2598), (15, 14, 0.1916), 
    (17, 14, 0.0867), (12, 15, 0.5706), (17, 15, 0.0267), (15, 16, 0.0859), 
    (8, 17, 0.011), (5, 18, 0.4737), (13, 18, 0.4188), (17, 18, 0.0526), 
    (7, 19, 0.66), (8, 19, 0.4883)
])
expected_max_real_2 = 0.4892312181238613

# =====================================================================
# EXAMPLE 3 DATA (Inverter PLL)
# =====================================================================
eqs_3 = [
    '-xp1 + (x3+ 0.2732*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1))+2*pi*50)*(-y2)',
    '-xp2 + (x3+ 0.2732*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1))+2*pi*50)*y1',
    '-xp3 + (0.022508948*1/3*((-2*u1+u2+u3)*x2+sqrt(3)*(u2-u3)*(x1)))',   
    'y1 - x1',
    'y2 - x2'
]
v_dict_3 = {
    'xp1': 0.0, 'xp2': 0.0, 'xp3': 0.0,
    'x1': -0.5737, 'x2': -0.8193, 'x3': -0.1597,
    'u1': -100.6109, 'u2': -217.5716, 'u3': 318.1824,
    'y1': -0.5737, 'y2': -0.8193
}

expected_E_3 = build_matrix((5, 5), [
    (0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)
])

expected_A_3 = build_matrix((5, 5), [
    (0, 0, -69.23543811590596), (1, 0, 48.4808627451425), (2, 0, -6.962416387610758), 
    (3, 0, -1.0), (0, 1, 22.520007971991998), (1, 1, -15.769228089261333), 
    (2, 1, 2.264644766034933), (4, 1, -1.0), (0, 2, 0.8193), (1, 2, -0.5737), 
    (1, 3, 339.9604201321298), (3, 3, 1.0), (0, 4, -339.9604201321298), (4, 4, 1.0)
])
expected_max_real_3 = -0.0058449811493328984


# =====================================================================
# FILLER DATA TO COMPLETE THE 5 EXAMPLES
# =====================================================================
eqs_4 = ["dx1 - x1 + u1"]
ineqs_4 = ["x1 - 10"]
v_dict_4 = {"dx1": 0.0, "x1": 5.0, "u1": 2.0}

eqs_5 = ["dx1 + 2*x1"]
v_dict_5 = {"dx1": 0.0, "x1": 1.0}


# =====================================================================
# TEST LIST DEFINITION
# Format: (Name, eqs, ineqs, v_dict, use_ineqs, expected_A, expected_E, expected_max_real)
# =====================================================================
example_list = [
    ("Example 1 - Fractions", eqs_1, [], v_dict_1, False, expected_A_1, expected_E_1, expected_max_real_1),
    ("Example 2 - 20-variable system", eqs_2, [], v_dict_2, False, expected_A_2, expected_E_2, expected_max_real_2),
    ("Example 3 - Inverter PLL", eqs_3, [], v_dict_3, False, expected_A_3, expected_E_3, expected_max_real_3),
    # Examples 4 and 5 run the base flow without strict matrix validation
    ("Example 4 - Inequalities", eqs_4, ineqs_4, v_dict_4, True, None, None, None), 
    ("Example 5 - Stable Linear System", eqs_5, [], v_dict_5, False, None, None, None),
]

# =====================================================================
# THE PARAMETRIZED TEST
# =====================================================================
@pytest.mark.parametrize(
    "name, eqs, ineqs, v_dict, use_ineqs, expected_A, expected_E, expected_max_real", 
    example_list
)
def test_builder_checks_matrices_and_stability(
    name, eqs, ineqs, v_dict, use_ineqs, expected_A, expected_E, expected_max_real
):
    # 1. Arrange & Act
    builder = PolynomialMatrixBuilder(eqs, ineqs, verbose=False)
    builder.linearize(v_dict, use_ineqs=use_ineqs)
    
    # 2. Assert: Check Matrix A 
    if expected_A is not None:
        np.testing.assert_allclose(
            builder.A, 
            expected_A, 
            atol=1e-4, 
            err_msg=f"{name}: Calculated matrix A does not match the expected one."
        )

    # 3. Assert: Check Matrix E 
    if expected_E is not None:
        np.testing.assert_allclose(
            builder.E, 
            expected_E, 
            atol=1e-4,
            err_msg=f"{name}: Calculated matrix E does not match the expected one."
        )

    # 4. Assert: Check Eigenvalue (max_real)
    if not use_ineqs:
        _, _, _, _, _, max_real = builder.compute_stability()
        
        if expected_max_real is None:
            assert max_real is None, f"{name}: Expected max_real=None, but got {max_real}"
        else:
            assert max_real == pytest.approx(expected_max_real, abs=1e-4), \
                f"{name}: Calculated max_real differs from the expected value."