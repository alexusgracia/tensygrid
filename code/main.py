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