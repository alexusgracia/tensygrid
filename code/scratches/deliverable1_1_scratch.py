import sympy as sp
import numpy as np
from scipy.linalg import block_diag
from sympy.parsing.sympy_parser import parse_expr


class PolynomialMatrixBuilder:
    """
    This class is used to build the S_H and Phi_F matrices from a list of polynomial equations.
    """

    def __init__(self, eqs, ineqs):
        """
        Initializes the PolynomialMatrixBuilder.
        
        Args:
            eqs (list): List of polynomial equations.
            ineqs (list): List of polynomial inequalities.
        """
        self.eqs = [parse_expr(eq, evaluate=False) for eq in eqs]
        self.ineqs = [parse_expr(eq, evaluate=False) for eq in ineqs]
        self.all_symbols = None

    def extract_symbols(self, all_exprs):
        """
        Extracts the symbols from a list of equations.
        
        Args:
            all_exprs (list): List of equations.
        
        Returns:
            list: List of symbols (sorted).
        """  
        all_symbols = []
        seen = set()
        for eq in all_exprs:
            # This is performed in order to have a consistent order of symbols
            current_syms = sorted(list(eq.free_symbols), key=lambda s: s.name)
            for sym in current_syms:
                if sym not in seen:
                    seen.add(sym)
                    all_symbols.append(sym)

        def custom_sort_key(sym):
            name = sym.name
            # This only orders the symbols as TenSyGrid deliverable 1.1 expects them
            if name.startswith('dx'): return (0, name)
            if name.startswith('x'): return (1, name)
            if name.startswith('u'): return (2, name)
            if name.startswith('y'): return (3, name)
            if name.startswith('z'): return (4, name)
            return (5, name)

        all_symbols.sort(key=custom_sort_key)

        print(f"Unified symbols order: {all_symbols}\n")
        self.all_symbols = all_symbols
        return all_symbols


    def matrix_creation(self, all_exprs):
        """
        Creates the S_H and Phi_F matrices from a list of equations.
        
        Args:
            all_exprs (list): List of equations.
            all_symbols (list): List of symbols.
        
        Returns:
            tuple: Tuple containing S_H and Phi_F matrices.
        """
        S_H_list = []
        P_Hs_list = []
        print(f"--- Expressions: ---")
        for i, eq in enumerate(all_exprs):
            print(f"\nEquation {i+1} (Ordered Terms): {eq}")
            try:
                if eq.is_Add:
                    terms = eq.args
                else:
                    terms = [eq]
                
                eq_coeffs = []
                
                for term in terms:
                    p = sp.Poly(term, *self.all_symbols)
                    monom = np.array(p.monoms()[0])
                    coeff = p.coeffs()[0]
                    
                    # Comment this line if you don't want to see the term comparison
                    print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom}")
                    
                    S_H_list.append(monom)
                    eq_coeffs.append(coeff)

                if eq_coeffs:
                    P_Hs_list.append(np.array(eq_coeffs))


            except sp.PolificationFailed:
                 print(f"\nCould not create Poly for equation {i+1}")
                
        if S_H_list:
            S_H = np.vstack(S_H_list).T
        else:
            S_H = np.array([])
            
        if P_Hs_list:
            # \Phi_F is a block diagonal matrix with the coefficients of each equation on the diagonal
            Phi_F = block_diag(*P_Hs_list)
        else:
            Phi_F = np.array([])
        
        return S_H, Phi_F


    def build(self):
        print("\n--- Coefficient extraction ---\n")
        self.extract_symbols(self.eqs + self.ineqs)

        S_H, Phi_F = self.matrix_creation(self.eqs)
        S_W, Phi_W = self.matrix_creation(self.ineqs)

        print("\n--- Matrix creation ---\n")
        print("S_H=\n", S_H)
        print("\n" + "\\Phi_F=" + "\n", Phi_F)
        print("\n" + "S_W=" + "\n", S_W)
        print("\n" + "\\Phi_W=" + "\n", Phi_W)

        return S_H, Phi_F, S_W, Phi_W



# --- Unfactorized case ---
### NOTE: This case works properly

eqs = [
    "3*dx1*y1 + 2*z1*x1 - u1*x1 - u1*z1*x1 + 2*x1",
    "dx1-y1"
]

ineqs = [
    "5 - 5*z1 - x1 + x1*z1"
]

builder = PolynomialMatrixBuilder(eqs, ineqs)
S_H, Phi_F, S_W, Phi_W = builder.build()


# --- Factorized case ---
### NOTE: This case does not work properly

eqs = [
    "3*dx1*y1 + 6*(1/2+1/2*z1)*(2/3-1/3*u1)*x1",
    "dx1-y1"
]

ineqs = [
    "12*(5/6-1/6*x1)*(1/2+1/2*z1)"
]

builder = PolynomialMatrixBuilder(eqs, ineqs)
S_H, Phi_F, S_W, Phi_W = builder.build()
