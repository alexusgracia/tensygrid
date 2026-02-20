import sympy as sp
import numpy as np
from scipy import sparse, linalg
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
            verbose (boolean): Returns debugging information.

        Returns:
            None
        """
        self.verbose = verbose
        self.eqs = [parse_expr(eq, evaluate=False) for eq in eqs]
        self.ineqs = [parse_expr(eq, evaluate=False) for eq in ineqs]
        self.all_symbols = None
        
        self.extract_symbols(self.eqs + self.ineqs)

        # Mapa de les variables simbòliques a nombre
        self.sym_to_idx = {sym.name: i for i, sym in enumerate(self.all_symbols)}
        
        S_H_sp, Phi_H_sp = self.matrix_creation(self.eqs)
        S_W_sp, Phi_W_sp = self.matrix_creation(self.ineqs)
        
        self.S_H = S_H_sp.toarray()
        self.Phi_H = Phi_H_sp.toarray()
        self.S_W = S_W_sp.toarray()
        self.Phi_W = Phi_W_sp.toarray()

        self.E = None
        self.A = None
        self.B = None

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
        Creates S_H and Phi_F matrices with automatic L1 normalization.

        Args:
            all_exprs (list): List of equations.
                
        Returns:
            S_H (np.array): Internal weights of variables (normalized b in a + b*s).
            Phi_F (np.array): Global coefficients of the terms (scaled by normalization factors).
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

        if S_H_list:
            S_H_raw = np.vstack(S_H_list).T.astype(float)
            self.S_H = sparse.csc_matrix(S_H_raw)
            
            num_cols = S_H_raw.shape[1]
            phi_rows = []
            for eq_dict in P_Hs_data:
                row = np.zeros(num_cols)
                for col_idx, val in eq_dict.items():
                    row[col_idx] = val
                phi_rows.append(row)
            
            self.Phi_H = sparse.csr_matrix(np.vstack(phi_rows))
        else:
            self.S_H = sparse.csc_matrix((len(self.all_symbols), 0))
            self.Phi_H = sparse.csr_matrix((len(all_exprs), 0))
        
        return self.S_H, self.Phi_H

    def linearize(self, v_dict, use_ineqs=False):
        """
        Analytic linearization for iMTI models (CPN representation).
        Matches the logic of the MATLAB scalable linearization.
        """
        S = self.S_W if use_ineqs else self.S_H
        Phi = self.Phi_W if use_ineqs else self.Phi_H
        
        v = np.zeros(len(self.all_symbols))
        for name, val in v_dict.items():
            if name in self.sym_to_idx:
                v[self.sym_to_idx[name]] = val

        X = (S.T * v) + (1 - np.abs(S.T)) 
        Y = np.prod(X, axis=1) 
        
        with np.errstate(divide='ignore', invalid='ignore'):
            invX = np.where(np.abs(X) > 1e-12, 1.0 / X, 0.0)
            F = S * (Y[:, np.newaxis] * invX).T
            
            zeros_per_column = np.sum(np.abs(X) < 1e-12, axis=1)
            for col_idx in np.where(zeros_per_column == 1)[0]:
                zero_row = np.where(np.abs(X[col_idx, :]) < 1e-12)[0][0]
                other_factors = np.delete(X[col_idx, :], zero_row)
                F[zero_row, col_idx] = S[zero_row, col_idx] * np.prod(other_factors)

        EABC = Phi @ F.T
        

        idx_dx = [i for i, s in enumerate(self.all_symbols) if s.name.startswith('dx')]
        idx_vars = [i for i, s in enumerate(self.all_symbols) if s.name.startswith(('x', 'y', 'z'))]
        idx_u = [i for i, s in enumerate(self.all_symbols) if s.name.startswith('u')]

        n_eqs = EABC.shape[0]
        n_vars = len(idx_vars)

        self.E = np.zeros((n_eqs, n_vars))
        for i, col_idx in enumerate(idx_dx):
            if i < n_vars:
                self.E[:, i] = -EABC[:, col_idx]

        self.A = EABC[:, idx_vars]
        self.B = EABC[:, idx_u]

        return EABC

    def compute_stability(self):
        """
        Computes the generalized eigenvalues for the pair (A, E).
        """
        if self.A is None or self.E is None:
            return None, False, None
        try:
            eigenvalues = linalg.eigvals(self.A, self.E)
            
            finite_evs = eigenvalues[np.isfinite(eigenvalues)]
            
            if len(finite_evs) == 0:
                return eigenvalues, True, -np.inf
            
            max_real = np.max(np.real(finite_evs))
            is_stable = max_real < 0
            
            return eigenvalues, is_stable, max_real
        except Exception as e:
            print(f"Error computing eigenvalues: {e}")
            return None, False, None

def profile_it(builder_instance, v_dict):
    """
    Profiles the execution of the linearize method.
    """
    pr = cProfile.Profile()
    pr.enable()
    
    builder_instance.linearize(v_dict)
    
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
    print("\n" + "="*40 + "\nTOP EXECUTION TIME\n" + "="*40)
    ps.print_stats(15)
    print(s.getvalue())

if __name__ == "__main__":
    do_profile = False 

    eqs = [
    "3*dx1*y1 + 6*(1/2+1/2*z1)*(2/3-1/3*u1)*x1",
    "y1 - 2*x1",
    "z1 + 0.5*y1 - 1"
]

    # ineqs =[
    #     "12*(5/6-1/6*x1)*(1/2-1/2*z1)"
    # ]
    builder = PolynomialMatrixBuilder(eqs, ineqs=[], verbose=False)
    
    v_dict = {
        'dx1': 1.0, 'x1': 2.0, 'u1': 3.0, 'y1': 4.0, 'z1': 5.0
    }
    

    print("\n--- S Matrix ---")
    print(builder.S_H)

    print("\n--- Phi Matrix ---")
    print(builder.Phi_H)

    if do_profile:
        profile_it(builder, v_dict)
    else:
        EABC = builder.linearize(v_dict)
        print("\n--- EABC ---")
        print(EABC)

    print("\n--- E Matrix ---")
    print(builder.E)

    evs, stable, margin = builder.compute_stability()
    print(f"\n--- Stability Analysis ---")
    print(f"Is stable: {stable}")
    print(f"Max Real Part: {margin}")