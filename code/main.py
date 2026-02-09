from MatrixBuilder import MatrixBuilder
import sympy as sp

equations_sets = [
    ['-dx1 - x1 + u','-dx2 - x2 + x1*x2',],
    ['-dx1 - x1 + u + 2*x2','-dx2 - x2 + x1*x2',],
    ['-3*dx1 + 2*x1 - 5*u','-dx2 - x2 + x1*x2',]
]

matrix = MatrixBuilder(equations_sets[1])
matrix.print_equations()

#print("S =")
#sp.pprint(builder.S) 
#print("P =")
#sp.pprint(builder.P)


### FORMATTING PROPOSAL ###
from sympy.parsing.latex import parse_latex
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
expr = r"-\frac{dx}{dt} - x1 + u+ 2x2"
transformations = (standard_transformations + (implicit_multiplication_application,))

equations = [
    '-dx1 - x1 + u+ 2x2',
    parse_latex(expr),
    sp.parse_expr("-dx1 - x1 + u+ 2x2", transformations=transformations),
    '-dx2 - x2 + x1*x2',

]

matrix = MatrixBuilder(equations)
matrix.print_equations()

