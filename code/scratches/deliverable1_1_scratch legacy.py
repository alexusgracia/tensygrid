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
        Creates the S_H and Phi_F matrices from a list of equations.
        
        Args:
            all_exprs (list): List of equations.
            all_symbols (list): List of symbols.
        
        Returns:
            tuple: Tuple containing S_H and Phi_F matrices.
        """
        S_H_list = []
        P_Hs_data = []
        monom_to_idx = {}
        
        if self.verbose:
            print(f"--- Expressions: ---")
            
        for i, eq in enumerate(all_exprs):
            if self.verbose:
                print(f"\nEquation {i+1} (Ordered Terms): {eq}")
            try:
                if eq.is_Add:
                    terms = eq.args
                else:
                    terms = [eq]
                
                current_eq_coeffs = {}
                
                for term in terms:
                    p = sp.Poly(term, *self.all_symbols)
                    monom_tuple = p.monoms()[0]
                    coeff = float(p.coeffs()[0])
                    
                    if monom_tuple not in monom_to_idx:
                        idx = len(S_H_list)
                        monom_to_idx[monom_tuple] = idx
                        S_H_list.append(np.array(monom_tuple))
                        
                        if self.verbose:
                            # Comment this line if you don't want to see the term comparison
                            print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom_tuple}")
                    else:
                        idx = monom_to_idx[monom_tuple]
                        if self.verbose:
                            print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom_tuple} already exists. Index = {idx}")
                    
                    current_eq_coeffs[idx] = current_eq_coeffs.get(idx, 0.0) + coeff

                P_Hs_data.append(current_eq_coeffs)

            except sp.PolificationFailed:
                if self.verbose:
                    print(f"\nCould not create Poly for equation {i+1}")
                
        if S_H_list:
            S_H = np.vstack(S_H_list).T
            num_cols = S_H.shape[1]
            
            phi_rows = []
            for eq_dict in P_Hs_data:
                row = np.zeros(num_cols)
                for col_idx, val in eq_dict.items():
                    row[col_idx] = val
                phi_rows.append(row)
            
            Phi_F = sparse.csr_matrix(np.vstack(phi_rows))
            S_H = sparse.csc_matrix(S_H)
        else:
            S_H = sparse.csc_matrix([])
            Phi_F = sparse.csc_matrix([])
        
        return S_H, Phi_F

    """
    @profile
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

        return EABC
        """

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
    do_profile = False  # <--- Cambia esto a True cuando quieras el reporte

    eqs = [
        "3*dx1*y1 + 2*z1*x1-u1*x1-u1*z1*x1+2*x1",
        "dx1-y1-x1"
    ]
    ineqs = [
        "5-5*z1-x1+x1*z1"
    ]

    builder = PolynomialMatrixBuilder(eqs, ineqs, verbose=True)

    v_dict = {
        'dx1': 1.0,
        'x1': 2.0,
        'u1': 3.0,
        'y1': 4.0,
        'z1': 5.0
    }
    
    if do_profile:
        profile_it(builder, v_dict)
    """
    else:
        result = builder.linearize(v_dict)
        print("\n--- Resultado de Linearización (EABC) ---")
        print(result)
    """
    