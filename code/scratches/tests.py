from deliverable1_1_scratch import PolynomialMatrixBuilder as PMB
import numpy as np

if __name__ == "__main__":
    eqs = [
        "3*dx1*y1 + 2*z1*x1-u1*x1-u1*z1*x1+2*x1",
        "dx1-y1-x1"
    ]
    ineqs = [
        "5-5*z1-x1+x1*z1"
    ]

    v_dict = {
        'dx1': 1.0,
        'x1': 2.0,
        'u1': 3.0,
        'y1': 4.0,
        'z1': 5.0
    }

    # Inicializar con verbose=True para ver el proceso en consola
    builder = PMB(eqs, ineqs, verbose=True)

    
    print("\n" + "="*30)
    print("S_H")
    print("="*30)
    print(builder.S_H)
    
    print("\n" + "="*30)
    print("Phi_F")
    print("="*30)
    print(builder.Phi_F)

"""
    result = builder.linearize(v_dict)
    print("\n--- EABC ---")
    print(result)
"""