import sympy as sp

# Definim els símbols (ara que ja sabem que podem tenir-los de VeraGrid)
#dx1, y1, z1, x1, u1 = sp.symbols('dx1 y1 z1 x1 u1')

# Definim les equacions que volem parsejar
equations = [
    '3*dx1*y1 + 2*z1*x1 - u1*x1 - u1*z1*x1 + 2*x1',
    'dx1 - y1',
    '5 - 5*z1 - x1 + x1*z1'
]

for equation in equations:
    # Extraiem els símbols de l'equació
    expr = sp.sympify(equation)

    # Separem els termes que volem
    terms = sp.Add.make_args(expr)
    
    # Per cada terme, separem els coeficients i les variables
    for term in terms:
        variables = term.atoms(sp.Symbol)
        print(variables)

        if len(variables) > 1:
            print(f" -> Combinació (Monomi): {term} | Variables: {variables}")
        elif len(variables) == 1:
            print(f" -> Element bàsic(escalar): {term} | Variable: {variables}")
        else:
            print(f" -> Constant: {term}")