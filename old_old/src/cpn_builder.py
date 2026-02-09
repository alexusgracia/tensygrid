import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from typing import List

class CPNBuilder:
    """
    Builds CPN representation matrices S and P from equations.
    Generic version: Does not distinguish between derivatives, states, or inputs.
    """
    
    def __init__(self, equations: List[str]):
        self.equations = equations
        self._parse_equations()
        self._order_symbols()       # Ahora es mucho más simple
        self._identify_monomials()
        self._build_matrices()
    
    def _parse_equations(self):
        """Parse equations and create a master dictionary of symbols."""
        # 1. Extract all symbol names generically
        all_symbol_names = set()
        for eq in self.equations:
            # Parse temporal para extraer nombres
            all_symbol_names.update([s.name for s in parse_expr(eq).free_symbols])
        
        # 2. Master dictionary (CRUCIAL for memory identity)
        # Esto asegura que el objeto Symbol('x') sea el mismo en todas partes
        self.master_subs = {name: sp.Symbol(name) for name in sorted(all_symbol_names)}
        
        # 3. Parsear ecuaciones reales usando el diccionario maestro
        self.expressions = [parse_expr(eq, local_dict=self.master_subs) for eq in self.equations]
    
    def _order_symbols(self):
        """
        Order symbols strictly alphabetically.
        Simplification: We don't care if it's dx, x, or u.
        """
        # Obtenemos los valores del diccionario ordenados por su nombre
        self.ordered_symbols = sorted(self.master_subs.values(), key=lambda s: s.name)
    
    def _identify_monomials(self):
        """Identify monomials ensuring linear ones match the symbol order."""
        all_monomials = set()
        for e in self.expressions:
            all_monomials.update(e.expand().as_coefficients_dict().keys())
        
        # 1. Monomios Lineales:
        # Iteramos sobre los símbolos ordenados y vemos si existen como monomio.
        # Esto garantiza que la matriz S tenga una diagonal de 1s para la parte lineal.
        linear_m = [sym for sym in self.ordered_symbols if sym in all_monomials]
        
        # 2. Monomios No Lineales (Productos):
        # Los detectamos por tener >1 símbolo libre y los ordenamos por texto
        nonlinear_m = sorted([m for m in all_monomials if len(m.free_symbols) > 1], key=str)
        
        self.ordered_monomials = linear_m + nonlinear_m
    
    def _build_matrices(self):
        """Build S and P using exact object comparison."""
        # Matrix S (Structure)
        S_data = []
        for sym in self.ordered_symbols:
            row = []
            for m in self.ordered_monomials:
                # La magia: 'sym in m.free_symbols' funciona perfecto gracias a master_subs
                if sym in m.free_symbols:
                    row.append(1)
                else:
                    row.append(0)
            S_data.append(row)
        self.S = sp.Matrix(S_data)
        
        # Matrix P (Coefficients)
        P_data = []
        for e in self.expressions:
            coeffs = e.expand().as_coefficients_dict()
            row = [coeffs.get(m, 0) for m in self.ordered_monomials]
            P_data.append(row)
        self.P = sp.Matrix(P_data)

    def get_S(self):
        """
        Return the structure matrix S.
        Rows = symbols, Columns = monomials; S[i,j] = 1 if symbol i appears in monomial j.
        """
        return self.S

    def get_P(self):
        """
        Return the parameter/coefficient matrix P.
        Rows = equations, Columns = monomials; P[i,j] = coefficient of monomial j in equation i.
        """
        return self.P

    def print_info(self):
        print("-" * 30)
        
        
        print("-" * 30)
        print("Matrix S (Structure):")
        print(f"Symbols (Rows S): {[s.name for s in self.ordered_symbols]}")
        # Print the row name (symbol) at the beginning of each S matrix row, preserving matrix format (no column names)
        row_names = [s.name for s in self.ordered_symbols]
        s_matrix = self.S.tolist()
        maxlen = max(len(name) for name in row_names)
        for name, row in zip(row_names, s_matrix):
            print(f"{name.rjust(maxlen)}: [ " + "  ".join(str(x) for x in row) + " ]")
        print("-" * 30)
        print("Matrix P (Coefficients):")
        print(f"Monomials (Cols S/P): {[str(m) for m in self.ordered_monomials]}")
        # Print P matrix with column (monomial) names, no row names
        col_names = [str(m) for m in self.ordered_monomials]
        p_matrix = self.P.tolist()
        maxlen = max(len(name) for name in col_names)
        header = "      " + "  ".join(name.rjust(maxlen) for name in col_names)
        print(header)
        for row in p_matrix:
            print("    [ " + "  ".join(str(x).rjust(maxlen) for x in row) + " ]")


