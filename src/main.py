from cpn_builder import CPNBuilder

# Equations (strings)
equations1 = [
    '2*x + y - x*u',
    'x - 3*y + z - u*z',
    '4*z - u',
]

equations2 = [
    '-dx1 - x1 + u',
    '-dx2 - x2 + x1*x2',
]

equations = equations2

# Build CPN representation
builder = CPNBuilder(equations)


builder.print_info()

# Extract matrices S and P
S = builder.get_S()
P = builder.get_P()