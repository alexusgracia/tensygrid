import sympy as sp
import numpy as np
from scipy import sparse
from sympy.parsing.sympy_parser import parse_expr
import cProfile, pstats, io
from memory_profiler import profile
import builtins

# Si 'profile' no existe en el entorno global, creamos uno que no haga nada
if 'profile' not in builtins.__dict__:
    def profile(func): 
        return func

class PolynomialMatrixBuilder:
    """
    This class is used to build the S_H and Phi_F matrices from a list of polynomial equations.
    """

    def __init__(self, eqs, ineqs, verbose=False):
        """
        Initializes the PolynomialMatrixBuilder.
        
        Args:
            eqs (list): List of polynomial equations.
            ineqs (list): List of polynomial inequalities.
        """
        self.verbose = verbose
        self.eqs = [parse_expr(eq, evaluate=False) for eq in eqs]
        self.ineqs = [parse_expr(eq, evaluate=False) for eq in ineqs]
        self.all_symbols = None
        
        self.extract_symbols(self.eqs + self.ineqs)

        # Mapa de les variables simbòliques a nombre
        self.sym_to_idx = {sym.name: i for i, sym in enumerate(self.all_symbols)}
        
        S_H_sp, Phi_F_sp = self.matrix_creation(self.eqs)
        S_W_sp, Phi_W_sp = self.matrix_creation(self.ineqs)
        
        self.S_H = S_H_sp.toarray()
        self.Phi_F = Phi_F_sp.toarray()
        self.S_W = S_W_sp.toarray()
        self.Phi_W = Phi_W_sp.toarray()

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

        if self.verbose:
            print(f"Unified symbols order: {all_symbols}\n")
        self.all_symbols = all_symbols
        return all_symbols

    def matrix_creation(self, all_exprs):
        """
        Creates S_H and Phi_F matrices with automatic L1 normalization (Eq. 4.18 - 4.21).
        S_H: Internal weights of variables (normalized b in a + b*s).
        Phi_F: Global coefficients of the terms (scaled by normalization factors).
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
                global_phi = float(coeff_sym)

                weights = [0.0] * len(self.all_symbols)
                
                # 2. Normalize each factor following TenSyGrid papers (Eq. 4.19)
                for f in factors:
                    f_vars = list(f.free_symbols)
                    if len(f_vars) == 1:
                        s = f_vars[0]
                        idx = self.sym_to_idx[s.name]
                        
                        # Extract internal coefficients b (slope) and a (intercept)
                        b_val = float(sp.diff(f, s))
                        a_val = float(f.subs(s, 0))
                        
                        # Apply L1 Normalization: Scale = |a| + |b|
                        scale = abs(a_val) + abs(b_val)
                        if scale == 0: scale = 1.0
                        
                        # Normalized weight goes to S, scale moves to Phi
                        normalized_b = b_val / scale
                        global_phi *= scale
                        weights[idx] = normalized_b
                        
                    elif len(f_vars) > 1:
                        # For expanded terms or complex factors, we use unit weights
                        for s in f_vars:
                            weights[self.sym_to_idx[s.name]] = 1.0
                    else:
                        # Pure numerical factors
                        global_phi *= float(f)

                monom_tuple = tuple(weights)
                
                # 3. Handle indexing and verbose logging
                if monom_tuple not in monom_to_idx:
                    idx = len(S_H_list)
                    monom_to_idx[monom_tuple] = idx
                    S_H_list.append(np.array(monom_tuple))
                    if self.verbose:
                        print(f"  Unseen term: {term}")
                        print(f"    -> Coeff (Phi): {global_phi}, Weights (S): {monom_tuple}")
                else:
                    idx = monom_to_idx[monom_tuple]
                    if self.verbose:
                        print(f"  Seen term: {term}")
                        print(f"    -> Adding {global_phi} to index {idx}")
                
                current_eq_coeffs[idx] = current_eq_coeffs.get(idx, 0.0) + global_phi

            P_Hs_data.append(current_eq_coeffs)

        # Construcción final de las matrices sparse
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

    """@profile
    def linearize(self, v_dict, use_ineqs=False):
        S = self.S_W if use_ineqs else self.S_H
        Phi = self.Phi_W if use_ineqs else self.Phi_F

        v = np.zeros((len(self.all_symbols), 1))
        for name, val in v_dict.items():
            if name in self.sym_to_idx:
                v[self.sym_to_idx[name]] = val

        X = S * v - np.abs(S)
        
        icrit = (X == -1)
        if np.any(icrit):
            raise ValueError("Specials to be implemented")
        
        X = np.where(S != 0, X + 1, 1.0)

        # axis=0 suma/multiplica hacia abajo (por columnas)
        Y = np.prod(X, axis=0) 

        F = S * Y * (1.0 / X) 

        EABC = Phi @ F.T

        return EABC"""

def profile_it(builder_instance, v_dict):
    pr = cProfile.Profile()
    pr.enable()
    
    builder_instance.linearize(v_dict)
    
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
    print("\n" + "="*40 + "\nTOP TIEMPO DE EJECUCIÓN\n" + "="*40)
    ps.print_stats(15)
    print(s.getvalue())

if __name__ == "__main__":
    do_profile = False 

    eqs = [
        "3*dx1*y1 + 2*z1*x1 - u_1*x1-u1*z1*x1+2*x1",
        "dx1-y1"
    ]

    ineqs =[
        "(5-x1)*(1-z1)"
    ]
    builder = PolynomialMatrixBuilder(eqs, ineqs, verbose=True)
    
    print("\n--- S_H Matrix ---")
    print(builder.S_H)
    print("\n--- Phi_H Matrix ---")
    print(builder.Phi_F)

    print("\n--- S_W Matrix ---")
    print(builder.S_W)
    print("\n--- Phi_W Matrix ---")
    print(builder.Phi_W)

    v_dict = {
        'dx1': 1.0, 'x1': 2.0, 'u1': 3.0, 'y1': 4.0, 'z1': 5.0
    }
    
    """if do_profile:
        profile_it(builder, v_dict)
    else:
        result = builder.linearize(v_dict)
        print("\n--- Resultado de Linearización (EABC) ---")
        print(result)"""