import unittest
import sys
import os
import numpy as np
import sympy as sp

# Configuración de rutas
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import PolynomialMatrixBuilder_class as PMB

class TestPolynomialStability(unittest.TestCase):

    def test_unstable_20x20_system(self):
        """Test 1: Linear Unstable System (20 equations)"""
        eqs = [
            '-xp1 + 0.3440*x10 + 0.4474',
            '-xp2 + 0.0038*x2 + 0.6549*x10 + 0.7845*x13 + 0.0040*x14 + 0.2282*x19 + 0.9734',
            '-xp3 + 0.6267*u5 + 0.0651*x12 + 0.3473*x16 + 0.8588',
            '-xp4 + 0.0598*u5 + 0.2887*x17 + 0.0327*x18 + 0.8685',
            '-xp5 + 0.0831*u4 + 0.38984*x15 + 0.1381*x18 + 0.2472',
            '-xp6 + 0.4737*x8 + 0.0856*x16 + 0.2286*x18 + 1.6085',
            '-xp7 + 0.0838*u5 + 0.2301*x11 + 0.2002*x12 + 0.8884',
            '-xp8 + 0.0998*u5 + 0.66*x9 + 0.0515*x10 + 0.5883',
            '-xp9 + 0.6785*u2 + 0.0110*x7 + 0.4883*x9 + 0.2754*x12 + 0.6573',
            '-xp10 + 0.0276*x1 + 0.2263*x19 + 0.3483',
            '-xp11 + 0.0549*u2 + 0.2598*x3 + 0.2733*x12 + 0.0994*x16 + 0.1905*x20 + 1.9020',
            '-xp12 + 0.0732*u5 + 0.0466*u6 + 0.8391',
            '-xp13 + 0.1264*u1 + 0.5706*x5 + 0.0839*x13 + 0.2854*x20 + 1.4471',
            '-xp14 + 0.4188*x8 + 0.4638',
            '-xp15 + 0.7539*u1 + 0.4065*x15 + 0.5479',
            '-xp16 + 0.0105*u2 + 0.1236*x1 + 0.1916*x4 + 0.0859*x6 + 0.1613*x18 + 1.0436',
            '-xp17 + 0.2431*u3 + 0.0820',
            '-xp18 + 0.0578*x2 + 0.0867*x4 + 0.0267*x5 + 0.0526*x8 + 2.5279',
            '-xp19 + 0.2644*u5 + 0.4517*x12 + 1.1073',
            '-xp20 + 0.5101*u4 + 0.3917'
        ]
        
        builder = PMB.PolynomialMatrixBuilder(eqs, ineqs=[], verbose=False)
        xOP = [-0.1668, -0.459, 0.5462, -1.3023, -0.4115, -0.5401, 1.3665, 0.7068, 0.6076, 0.9198,
               0.8994, 0.0327, -1.4159, -1.6809, 0.0378, 1.6635, 0.831, 1.0477, 1.7938, -0.4851]
        uOP = [-0.6336, 1.1019, -0.7170, 2.2000, 0.8906, -0.6354]

        v_dict = {
            **{f'x{i+1}': xOP[i] for i in range(20)},
            **{f'u{i+1}': uOP[i] for i in range(len(uOP))},
            **{f'xp{i+1}': 0.0 for i in range(20)},
            **{f'dx{i+1}': 0.0 for i in range(20)}
        }

        builder.linearize(v_dict)
        _, _, _, _, stable, margin = builder.compute_stability()
        
        self.assertFalse(stable)
        self.assertAlmostEqual(margin, 0.489, places=3)

    def test_stable_nonlinear_5x5_system(self):
        """Test 2: Nonlinear Stable System (3 ODEs + 2 Algebraic)"""
        from math import pi
        sqrt = sp.sqrt
        
        # Note: Pi and sqrt(3) are parsed correctly by builder if string is clean
        eqs = [
            '-xp1 + (x3 + 0.2732*1/3*((-2*u1+u2+u3)*x2 + 1.73205081*(u2-u3)*x1) + 314.159265)*(-y2)',
            '-xp2 + (x3 + 0.2732*1/3*((-2*u1+u2+u3)*x2 + 1.73205081*(u2-u3)*x1) + 314.159265)*y1',
            '-xp3 + (0.022508948*1/3*((-2*u1+u2+u3)*x2 + 1.73205081*(u2-u3)*x1))',
            'y1 - x1',
            'y2 - x2'
        ]

        builder = PMB.PolynomialMatrixBuilder(eqs, ineqs=[], verbose=False)
        
        xOP = [-0.5737, -0.8193, -0.1597]
        uOP = [-100.6109, -217.5716, 318.1824]
        yOP = [-0.5737, -0.8193]

        v_dict = {
            **{f'x{i+1}': xOP[i] for i in range(len(xOP))},
            **{f'u{i+1}': uOP[i] for i in range(len(uOP))},
            **{f'y{i+1}': yOP[i] for i in range(len(yOP))},
            **{f'xp{i+1}': 0.0 for i in range(3)},
            **{f'dx{i+1}': 0.0 for i in range(3)}
        }

        builder.linearize(v_dict)
        evals, _, _, _, stable, margin = builder.compute_stability()

        self.assertTrue(stable)
        self.assertLess(margin, 0)
        self.assertAlmostEqual(margin, -0.00584, places=4)

if __name__ == "__main__":
    unittest.main()