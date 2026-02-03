from cpn_builder import CPNBuilder
from cpn_linearizer import CPNLinearizer

equations1 = [
    '-dx1 - x1 + u',
    '-dx2 - x2 + x1*x2',
]

equations2 = [
    '-dx1 - x1 + u + 2*x2',
    '-dx2 - x2 + x1*x2',
]

equations3 = [
    '-3*dx1 + 2*x1 - 5*u',
    '-dx2 - x2 + x1*x2',
]

equations = equations3

# Build CPN representation
builder = CPNBuilder(equations)


builder.print_info()


# Define parameters
n_states = 2
m_inputs = 1
p_outputs = 0
q_equations = n_states
N_signals = 2 * n_states + m_inputs + p_outputs

# Extract matrices S and P
S = builder.get_S()
P = builder.get_P()

linearizer = CPNLinearizer(equations, S, P, n_states, m_inputs, p_outputs, q_equations, N_signals)


linearizer.linearize(test_zeros = False, debug = False, total_time = False, print_eigenvalues = True)