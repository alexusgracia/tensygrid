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
                P_Hs_data = []
                monom_to_idx = {}
                
                print(f"--- Expressions: ---")
                for i, eq in enumerate(all_exprs):
                    print(f"\nEquation {i+1} (Ordered Terms): {eq}")
                    try:
                        if eq.is_Add:
                            terms = eq.args
                        else:
                            terms = [eq]
                        
                        # Diccionario local para esta ecuación: {indice_columna: suma_de_coeficientes}
                        current_eq_coeffs = {}
                        
                        for term in terms:
                            p = sp.Poly(term, *self.all_symbols)
                            monom_tuple = p.monoms()[0]
                            coeff = float(p.coeffs()[0])
                            
                            # Si el monomio no existe en la matriz S_H global, lo añadimos
                            if monom_tuple not in monom_to_idx:
                                idx = len(S_H_list)
                                monom_to_idx[monom_tuple] = idx
                                S_H_list.append(np.array(monom_tuple))
                                
                                # Comment this line if you don't want to see the term comparison
                                print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom_tuple}")
                            else:
                                idx = monom_to_idx[monom_tuple]
                                # Usamos tu comentario para mostrar que ya existe y su índice
                                print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom_tuple} already exists. Index = {idx}")
                            
                            # Acumulamos el coeficiente en la columna correspondiente para esta fila
                            current_eq_coeffs[idx] = current_eq_coeffs.get(idx, 0.0) + coeff

                        P_Hs_data.append(current_eq_coeffs)

                    except sp.PolificationFailed:
                        print(f"\nCould not create Poly for equation {i+1}")
                        
                if S_H_list:
                    # S_H se construye con los monomios únicos (uno por columna)
                    S_H = np.vstack(S_H_list).T
                    num_cols = S_H.shape[1]
                    
                    # Construimos Phi_F asegurando que cada fila tenga el ancho total de S_H
                    phi_rows = []
                    for eq_dict in P_Hs_data:
                        row = np.zeros(num_cols)
                        for col_idx, val in eq_dict.items():
                            row[col_idx] = val
                        phi_rows.append(row)
                    
                    # Convertimos a float para que scipy.sparse no falle
                    Phi_F = sparse.csr_matrix(np.vstack(phi_rows))
                    # Convertimos a sparse
                    S_H = sparse.csc_matrix(S_H)
                else:
                    S_H = sparse.csc_matrix([])
                    Phi_F = sparse.csc_matrix([])
                
                return S_H, Phi_F

    def linearize(self, S, Phi, v_dict):
        v = np.array([v_dict.get(sym, 0.0) for sym in self.all_symbols], dtype=float).reshape(-1, 1)

        

    @profile
    def build(self):
        print("\n--- Coefficient extraction ---\n")
        self.extract_symbols(self.eqs + self.ineqs)

        S_H, Phi_F = self.matrix_creation(self.eqs)
        S_W, Phi_W = self.matrix_creation(self.ineqs)
        S_H_array, Phi_F_array = S_H.toarray(), Phi_F.toarray()
        S_W_array, Phi_W_array = S_W.toarray(), Phi_W.toarray()
        print("\n--- Matrix creation ---\n")
        print("S_H=\n", S_H_array)
        print("\\Phi_F=\n", Phi_F_array)
        print("S_W=\n", S_W_array)
        print("\\Phi_W=\n", Phi_W_array)

        return S_H_array, Phi_F_array, S_W_array, Phi_W_array

# --- Orquestación de Profiling ---

def profile_it(builder_instance):
    # Tiempo
    pr = cProfile.Profile()
    pr.enable()
    
    # Esto disparará el @profile de memoria línea a línea en build()
    builder_instance.build()
    
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

    builder = PolynomialMatrixBuilder(eqs, ineqs)

    if do_profile:
        profile_it(builder)
    else:
        # Ejecución normal, limpia y rápida
        S_H, Phi_F, S_W, Phi_W = builder.build()
