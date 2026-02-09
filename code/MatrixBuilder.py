import sympy as sp

class MatrixBuilder:

    def __init__(self, equations=None):
        self.equations = equations
        self.S = None
        self.P = None
    
    def print_equations(self):
        for eq in self.equations:
            print(eq)