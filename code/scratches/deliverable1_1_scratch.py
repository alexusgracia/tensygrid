import sympy as sp
import numpy as np
from scipy import sparse
from sympy.parsing.sympy_parser import parse_expr
import cProfile, pstats, io
import builtins

# If 'profile' does not exist in the global environment, create a dummy decorator
if 'profile' not in builtins.__dict__:
    def profile(func): 
        return func

class PolynomialMatrixBuilder:
    """
    This class builds the S_H and Phi_F matrices from a list of polynomial equations
    following the TenSyGrid D1.1 Multilinear format.
    """

    def __init__(self, eqs, ineqs, verbose=False):
        """
        Initializes the builder and creates the matrices.
        
        Args:
            eqs (list): List of polynomial equations (strings).
            ineqs (list): List of polynomial inequalities (strings).
            verbose (bool): If True, prints the analysis of terms.
        """
        self.verbose = verbose
        # Parse expressions without evaluating to preserve user-defined factors
        self.eqs = [parse_expr(eq, evaluate=False) for eq in eqs]
        self.ineqs = [parse_expr(eq, evaluate=False) for eq in ineqs]
        self.all_symbols = None
        
        # Extract and sort symbols based on D1.1 naming conventions
        self.extract_symbols(self.eqs + self.ineqs)

        # Map symbolic variables to numeric indices
        self.sym_to_idx = {sym.name: i for i, sym in enumerate(self.all_symbols)}
        
        # Create matrices for equations and inequalities
        S_H_sp, Phi_F_sp = self.matrix_creation(self.eqs)
        S_W_sp, Phi_W_sp = self.matrix_creation(self.ineqs)
        
        # Convert sparse to dense for easier inspection
        self.S_H = S_H_sp.toarray()
        self.Phi_F = Phi_F_sp.toarray()
        self.S_W = S_W_sp.toarray()
        self.Phi_W = Phi_W_sp.toarray()

    def extract_symbols(self, all_exprs):
        """
        Extracts symbols from expressions and sorts them by category (dx, x, u, y, z).
        """  
        all_symbols = []
        seen = set()
        for eq in all_exprs:
            # Sort local symbols to maintain consistency
            current_syms = sorted(list(eq.free_symbols), key=lambda s: s.name)
            for sym in current_syms:
                if sym not in seen:
                    seen.add(sym)
                    all_symbols.append(sym)

        def custom_sort_key(sym):
            name = sym.name
            # Sorting logic as expected by TenSyGrid Deliverable 1.1
            if name.startswith('dx'): return (0, name)
            if name.startswith('x'): return (1, name)
            if name.startswith('u'): return (2, name)
            if name.startswith('y'): return (3, name)
            if name.startswith('z'): return (4, name)
            return (5, name)

        all_symbols.sort(key=custom_sort_key)

        if self.verbose:
            print(f"Unified symbols order: {all_symbols}\n")
        self.all_symbols = all_symbols
        return all_symbols

    def matrix_creation(self, all_exprs):
        """
        Creates S_H and Phi_F matrices.
        S_H: Internal weights of variables (b in a + b*s).
        Phi_F: Global coefficients of the terms.
        """
        S_H_list = []
        P_Hs_data = []
        monom_to_idx = {}
        
        if self.verbose:
            print(f"--- Analyzing Expressions for Matrices ---")
            
        for i, eq in enumerate(all_exprs):
            if self.verbose:
                print(f"\nEquation {i+1}: {eq}")
            
            # Treat the equation as a sum of terms
            terms = eq.args if eq.is_Add else [eq]
            
            current_eq_coeffs = {}
            for term in terms:
                # 1. Extract global coefficient and product factors
                coeff_sym, factors = term.as_coeff_mul()
                coeff = float(coeff_sym)

                # 2. Extract weights b_i for each symbol in this term
                weights = []
                for s in self.all_symbols:
                    s_weight = 0.0
                    if s in term.free_symbols:
                        # Search for the specific factor containing variable 's'
                        for f in factors:
                            if s in f.free_symbols:
                                # Differentiate to get internal coefficient 'b'
                                deriv = sp.diff(f, s)
                                # Clean residual variables for expanded cases (e.g., u1*z1)
                                for other in deriv.free_symbols:
                                    deriv = deriv.subs(other, 1)
                                s_weight = float(deriv)
                                break
                    weights.append(s_weight)
                
                monom_tuple = tuple(weights)
                
                # 3. Handle indexing and verbose logging
                if monom_tuple not in monom_to_idx:
                    idx = len(S_H_list)
                    monom_to_idx[monom_tuple] = idx
                    S_H_list.append(np.array(monom_tuple))
                    if self.verbose:
                        print(f"  Unseen term: {term}")
                        print(f"    -> Coeff (Phi): {coeff}, Weights (S): {monom_tuple}")
                else:
                    idx = monom_to_idx[monom_tuple]
                    if self.verbose:
                        print(f"  Seen term: {term}")
                        print(f"    -> Adding {coeff} to index {idx}")
                
                current_eq_coeffs[idx] = current_eq_coeffs.get(idx, 0.0) + coeff

            P_Hs_data.append(current_eq_coeffs)

        # Build final sparse matrices
        if S_H_list:
            S_H_raw = np.vstack(S_H_list).T.astype(float)
            S_H = sparse.csc_matrix(S_H_raw)
            
            num_cols = S_H_raw.shape[1]
            phi_rows = []
            for eq_dict in P_Hs_data:
                row = np.zeros(num_cols)
                for col_idx, val in eq_dict.items():
                    row[col_idx] = val
                phi_rows.append(row)
            
            Phi_F = sparse.csr_matrix(np.vstack(phi_rows))
        else:
            S_H = sparse.csc_matrix((len(self.all_symbols), 0))
            Phi_F = sparse.csr_matrix((len(all_exprs), 0))
        
        return S_H, Phi_F

if __name__ == "__main__":
    # Example equations combining factored and expanded forms
    eqs = [
        "3*dx1*y1 + 6*(1/2+1/2*z1)*(2/3-1/3*u1)*x1",
        "dx1-y1+dx1*y1"
    ]

    # Initialize builder with verbose=True to see the analysis
    builder = PolynomialMatrixBuilder(eqs, ineqs=[], verbose=True)
    
    print("\n--- Final S_H Matrix ---")
    print(builder.S_H)
    print("\n--- Final Phi_F Matrix ---")
    print(builder.Phi_F)