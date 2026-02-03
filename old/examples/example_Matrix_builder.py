from MatrixBuilder import MatrixBuilder
import sympy as sp

equations = [
    '-dx1 - x1 + u+ 2x2',
    '-dx2 - x2 + x1*x2',
]

builder = MatrixBuilder(equations)


print("S =")
sp.pprint(builder.S) 
print("P =")
sp.pprint(builder.P)


