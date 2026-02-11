import sympy as sp
import numpy as np
from scipy.linalg import block_diag
from sympy.parsing.sympy_parser import parse_expr


def extract_symbols(all_exprs):
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
        current_syms = sorted(list(eq.free_symbols), key=lambda s: s.name)
        for sym in current_syms:
            if sym not in seen:
                seen.add(sym)
                all_symbols.append(sym)

    def custom_sort_key(sym):
        name = sym.name
        if name.startswith('dx'): return (0, name)
        if name.startswith('x'): return (1, name)
        if name.startswith('u'): return (2, name)
        if name.startswith('y'): return (3, name)
        if name.startswith('z'): return (4, name)
        return (5, name)

    all_symbols.sort(key=custom_sort_key)

    print(f"Unified symbols order: {all_symbols}\n")
    return all_symbols



def matrix_creation(all_exprs, all_symbols):
    """
    Creates the S_H and Phi_F matrices from a list of equations.
    
    Args:
        all_exprs (list): List of equations.
        all_symbols (list): List of symbols.
    
    Returns:
        tuple: Tuple containing S_H and Phi_F matrices.
    """
    S_H_list = []
    P_Hs_list = []
    print(f"--- Expressions: ---")
    for i, eq in enumerate(all_exprs):
        print(f"\nEquation {i+1} (Ordered Terms): {eq}")
        try:
            if eq.is_Add:
                terms = eq.args
            else:
                terms = [eq]
            
            eq_coeffs = []
            
            for term in terms:
                p = sp.Poly(term, *all_symbols)
                monom = np.array(p.monoms()[0])
                coeff = p.coeffs()[0]
                
                print(f"  Term: {term} -> Coeff: {coeff}, Monom: {monom}")
                
                S_H_list.append(monom)
                eq_coeffs.append(coeff)

            if eq_coeffs:
                P_Hs_list.append(np.array(eq_coeffs))


        except sp.PolificationFailed:
             print(f"\nCould not create Poly for equation {i+1}")
            
    if S_H_list:
        S_H = np.vstack(S_H_list).T
    else:
        S_H = np.array([])
        
    if P_Hs_list:
        Phi_F = block_diag(*P_Hs_list)
    else:
        Phi_F = np.array([])
    
    return S_H, Phi_F

eq_1 = "3*dx1*y1 + 2*z1*x1 - u1*x1 - u1*z1*x1 + 2*x1"
eq_2 = "dx1-y1"
eq_3 = "5 - 5*z1 - x1 +x1*z1"

eq_1 = parse_expr(eq_1, evaluate=False)
eq_2 = parse_expr(eq_2, evaluate=False)
eq_3 = parse_expr(eq_3, evaluate=False)

eqs = [eq_1, eq_2]
ineqs = [eq_3]

print("\n--- Coefficient extraction ---\n")
all_symbols = extract_symbols(eqs + ineqs)
S_H, Phi_F = matrix_creation(eqs, all_symbols)
S_H_ineq, Phi_F_ineq = matrix_creation(ineqs, all_symbols)
print("\n--- Matrix creation ---\n")
print("S_H=\n", S_H)
print("\n" + "\\Phi_F=" + "\n", Phi_F)
print("\n" + "S_H_ineq=" + "\n", S_H_ineq)
print("\n" + "\\Phi_F_ineq=" + "\n", Phi_F_ineq)
