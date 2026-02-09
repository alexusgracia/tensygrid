from MatrixBuilder import MatrixBuilder
import sympy as sp

equations_sets = [
    ['-dx1 - x1 + u','-dx2 - x2 + x1*x2',],
    ['-dx1 - x1 + u + 2*x2','-dx2 - x2 + x1*x2',],
    ['-3*dx1 + 2*x1 - 5*u','-dx2 - x2 + x1*x2',],
    ['x1*(x2+x3*u1)*(x4+x5*y2)*(x6+x7*y3) + x8*(x9+x10*u1)*(x11+x12*y2)'],
    ['mw*cw*xp1 = cw*u1*(u2-x1)+Aw*u4*(u3-x1)', "0 = Aw*u1*(u2-y1) + Aa*(u4-y1)*u3"],

]

build = MatrixBuilder(equations_sets[3])
build.print_equations()

print("S =")
sp.pprint(build.S) 
print("P =")
sp.pprint(build.P)