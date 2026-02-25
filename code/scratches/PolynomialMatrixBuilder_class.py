import sympy as sp
import numpy as np
from scipy import sparse, linalg
from sympy.parsing.sympy_parser import parse_expr
from memory_profiler import profile
import builtins

if 'profile' not in builtins.__dict__:
    def profile(func):
        return func


class PolynomialMatrixBuilder:
    """
    Builds the S_H / Phi_H (and S_W / Phi_W) matrices from a list of
    polynomial equations, and performs analytic linearization using the
    iMTI / CPN representation described in the TenSyGrid papers.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, eqs: list[str], ineqs: list[str], verbose: bool = False):
        """
        Args:
            eqs     (list[str]): Polynomial equality constraints.
            ineqs   (list[str]): Polynomial inequality constraints.
            verbose (bool)     : Print debugging information.
        """
        self.verbose = verbose

        self.eqs   = [parse_expr(eq, evaluate=False) for eq in eqs]
        self.ineqs = [parse_expr(eq, evaluate=False) for eq in ineqs]

        self.all_symbols: list = []
        self.sym_to_idx:  dict = {}
        self.extract_symbols(self.eqs + self.ineqs)
        self.sym_to_idx = {sym.name: i for i, sym in enumerate(self.all_symbols)}

        S_H_sp,  Phi_H_sp  = self.matrix_creation(self.eqs)
        S_W_sp,  Phi_W_sp  = self.matrix_creation(self.ineqs)

        self.S_H   = S_H_sp.toarray()
        self.Phi_H = Phi_H_sp.toarray()
        self.S_W   = S_W_sp.toarray()
        self.Phi_W = Phi_W_sp.toarray()

        self.E: np.ndarray | None = None
        self.A: np.ndarray | None = None
        self.B: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Symbol extraction
    # ------------------------------------------------------------------

    def extract_symbols(self, all_exprs: list) -> list:
        """
        Extracts and sorts symbols from a list of SymPy expressions.

        Ordering convention (TenSyGrid deliverable 1.1):
            dx* < x* < u* < y* < z* < everything else

        Args:
            all_exprs (list): SymPy expressions to scan.

        Returns:
            list: Sorted list of unique symbols.
        """
        seen: set   = set()
        all_symbols = []
        for eq in all_exprs:
            for sym in sorted(eq.free_symbols, key=lambda s: s.name):
                if sym not in seen:
                    seen.add(sym)
                    all_symbols.append(sym)

        def _sort_key(sym):
            name = sym.name
            if name.startswith('dx'): return (0, name)
            if name.startswith('xp'): return (0, name)
            if name.startswith('x') and not name.startswith('xp'):  return (1, name)
            if name.startswith('u'):  return (2, name)
            if name.startswith('y'):  return (3, name)
            if name.startswith('z'):  return (4, name)
            return (5, name)

        all_symbols.sort(key=_sort_key)

        if self.verbose:
            print(f"Unified symbols order: {all_symbols}\n")

        self.all_symbols = all_symbols
        return all_symbols

    # ------------------------------------------------------------------
    # Matrix creation
    # ------------------------------------------------------------------

    @profile
    def matrix_creation(self, all_exprs: list) -> tuple[sparse.csc_matrix, sparse.csr_matrix]:
        """
        Creates the S and Phi sparse matrices with automatic L1 normalisation
        (TenSyGrid Eq. 4.19).

        Each column of S holds the normalised internal weights for one
        monomial term; each row of Phi holds the global coefficients of the
        monomials that appear in one equation.

        Args:
            all_exprs (list): SymPy expressions to encode.

        Returns:
            tuple: (S, Phi) - both sparse matrices.
        """
        S_list:       list = []
        Phi_data:     list = []
        monom_to_idx: dict = {}

        if self.verbose:
            print("--- Analyzing Expressions for Matrices ---")

        for i, eq in enumerate(all_exprs):
            if self.verbose:
                print(f"\nEquation {i + 1}: {eq}")

            terms = eq.args if eq.is_Add else [eq]
            current_eq_coeffs: dict = {}

            for term in terms:
                current_eq_coeffs = self._get_monomial_weights(term, S_list, monom_to_idx, current_eq_coeffs)
                

            Phi_data.append(current_eq_coeffs)

        # 4. Assemble sparse matrices
        if S_list:
            S_raw    = np.vstack(S_list).T.astype(float)
            S        = sparse.csc_matrix(S_raw)
            num_cols = S_raw.shape[1]
            phi_rows = []
            for eq_dict in Phi_data:
                row = np.zeros(num_cols)
                for col_idx, val in eq_dict.items():
                    row[col_idx] = val
                phi_rows.append(row)
            Phi = sparse.csr_matrix(np.vstack(phi_rows))
        else:
            S   = sparse.csc_matrix((len(self.all_symbols), 0))
            Phi = sparse.csr_matrix((len(all_exprs), 0))

        return S, Phi

    # ------------------------------------------------------------------
    # Linearisation – public entry point
    # ------------------------------------------------------------------

    @profile
    def linearize(
        self, v_dict: dict[str, float], use_ineqs: bool = False
    ) -> np.ndarray:
        """
        Analytic linearization for iMTI models (CPN representation).

        Computes the Jacobian of the polynomial system at the operating
        point defined by *v_dict*, then splits the result into the
        descriptor matrices E, A, B stored on the instance.

        Args:
            v_dict    (dict): Variable name → value at the operating point.
            use_ineqs (bool): Use the inequality matrices instead.

        Returns:
            np.ndarray: Full EABC Jacobian matrix.
        """
        S, Phi       = self._get_matrices(use_ineqs)
        v            = self._build_v_vector(v_dict)
        F            = self._compute_jacobian(S, v)
        EABC         = Phi @ F.T
        self.E, self.A, self.B = self._split_EABC(EABC)
        return EABC

    # ------------------------------------------------------------------
    # Linearisation – private helpers
    # ------------------------------------------------------------------

    def _get_matrices(self, use_ineqs: bool) -> tuple[np.ndarray, np.ndarray]:
        """Return the (S, Phi) pair for equations or inequalities."""
        if use_ineqs:
            return self.S_W, self.Phi_W
        return self.S_H, self.Phi_H

    def _build_v_vector(self, v_dict: dict[str, float]) -> np.ndarray:
        """Build the operating-point vector from a name→value dictionary."""
        v = np.zeros(len(self.all_symbols))
        for name, val in v_dict.items():
            if name in self.sym_to_idx:
                v[self.sym_to_idx[name]] = val
        return v

    def _get_monomial_weights(self, term, S_list, monom_to_idx, current_eq_coeffs):
            """
            Processes a term. If the term is non-linear/factored (e.g., (x+y)*z),
            it catches the error and expands it automatically.
            """
            try:
                # 1. Separate global coefficient
                coeff_sym, symbolic_part = term.as_coeff_mul()
                global_phi = float(coeff_sym)
                weights = [0.0] * len(self.all_symbols)

                # 2. Extract symbolic factors
                factors = sp.Mul.make_args(sp.Mul(*symbolic_part))

                for f in factors:
                    f_vars = list(f.free_symbols)

                    if len(f_vars) == 1:
                        s = f_vars[0]
                        idx = self.sym_to_idx[s.name]

                        # Extraction logic (b*s + a)
                        b_val = float(sp.diff(f, s))
                        a_val = float(f.subs(s, 0))

                        scale = abs(a_val) + abs(b_val)
                        if scale == 0: scale = 1.0

                        weights[idx] = b_val / scale
                        global_phi *= scale

                    elif len(f_vars) > 1:
                        # TRIGGER: This is a factored term like (x + y + 1)
                        # We raise a ValueError to trigger the 'except' block expansion
                        raise ValueError("Multi-variable factor detected.")
                    
                    elif len(f_vars) == 0:
                        global_phi *= float(f)

                # 3. Standard Registration (if no error was raised)
                monom_tuple = tuple(weights)
                if monom_tuple not in monom_to_idx:
                    idx = len(S_list)
                    monom_to_idx[monom_tuple] = idx
                    S_list.append(np.array(monom_tuple))
                else:
                    idx = monom_to_idx[monom_tuple]

                current_eq_coeffs[idx] = current_eq_coeffs.get(idx, 0.0) + global_phi

            except (ValueError, TypeError):
                # SELF-HEALING: Expand the problematic term and process its parts
                if self.verbose:
                    print(f"   Expanding complex term: {term}")
                
                expanded_sub_terms = sp.Add.make_args(term.expand())
                for sub_t in expanded_sub_terms:
                    # Recurse with simpler terms
                    self._get_monomial_weights(sub_t, S_list, monom_to_idx, current_eq_coeffs)

            return current_eq_coeffs

    @profile
    def _compute_jacobian(self, S: np.ndarray, v: np.ndarray) -> np.ndarray:
        """
        Compute the analytic Jacobian of the monomial basis at point *v*.

        For each monomial m(v) = ∏ᵢ (sᵢ·vᵢ + (1-|sᵢ|)), the partial
        derivative with respect to vⱼ is  sⱼ · ∏_{i≠j} Xᵢ  = sⱼ · Y / Xⱼ.

        Args:
            S (np.ndarray): Weight matrix (n_vars x n_monomials).
            v (np.ndarray): Operating-point vector (n_vars,).

        Returns:
            np.ndarray: Jacobian F (n_vars x n_monomials).
        """
        X = (S.T * v) + (1 - np.abs(S.T))   # shape: (n_monomials, n_vars)
        Y = np.prod(X, axis=1)               # shape: (n_monomials,)

        with np.errstate(divide='ignore', invalid='ignore'):
            invX = np.where(np.abs(X) > 1e-12, 1.0 / X, 0.0)
            F    = S * (Y[:, np.newaxis] * invX).T

            # Correct columns where exactly one factor is zero
            zeros_per_col = np.sum(np.abs(X) < 1e-12, axis=1)
            for col in np.where(zeros_per_col == 1)[0]:
                zero_row         = np.where(np.abs(X[col, :]) < 1e-12)[0][0]
                other_factors    = np.delete(X[col, :], zero_row)
                F[zero_row, col] = S[zero_row, col] * np.prod(other_factors)

        return F

    def _split_EABC(
        self, EABC: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Split the full Jacobian into the descriptor system matrices E, A, B.

        Args:
            EABC (np.ndarray): Full Jacobian (n_eqs x n_vars).

        Returns:
            tuple: (E, A, B) matrices.
        """
        idx_dx   = [i for i, s in enumerate(self.all_symbols)
                    if s.name.startswith('dx') or s.name.startswith('xp')]
        idx_vars = [i for i, s in enumerate(self.all_symbols)
                    if s.name.startswith(('x', 'y', 'z'))
                    and not s.name.startswith('xp')]
        idx_u    = [i for i, s in enumerate(self.all_symbols) if s.name.startswith('u')]

        n_vars = len(idx_vars)
        E      = np.zeros((EABC.shape[0], n_vars))
        for i, col in enumerate(idx_dx):
            if i < n_vars:
                E[:, i] = -EABC[:, col]

        A = EABC[:, idx_vars]
        B = EABC[:, idx_u]

        return E, A, B

    # ------------------------------------------------------------------
    # Stability analysis
    # ------------------------------------------------------------------

    @profile
    def compute_stability(
        self,
    ) -> tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None, np.ndarray | None, bool, float | None]:
        """
        Compute generalised eigenvalues for the descriptor pair (A, E).

        Returns:
            tuple: (eigenvalues, is_stable, max_real_part)
                   All entries are None / False on failure.
        """
        if self.A is None or self.E is None:
            return None, None, None, None, False, None

        try:
            evals, evecs_left, evecs_right = linalg.eig(self.A, self.E, left=True, right=True)
            finite_evs  = evals[np.isfinite(evals)]

            if len(finite_evs) == 0:
                return evals, evecs_left, evecs_right, None, True, -np.inf

            for i in range(len(evals)):
                scaling = np.dot(evecs_left[:, i].conj(), evecs_right[:, i])
                evecs_left[:, i] /= scaling

            participation_matrix = np.abs(evecs_right) * np.abs(evecs_left).T

            max_real  = float(np.max(np.real(finite_evs)))
            is_stable = max_real < 0
            return evals, evecs_left, evecs_right, participation_matrix, is_stable, max_real

        except Exception as e:
            print(f"Error computing eigenvalues: {e}")
            if self.verbose:
                try:
                    print(f"det(E - A) = {linalg.det(self.E - self.A)}")
                except Exception:
                    print("Error computing det(E - A)")
            return None, None, None, None, False, None

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def report(
        self,
        eigenvalues: np.ndarray | None = None,
        is_stable: bool = False,
        max_real: float | None = None,
        print_matrices: bool = True,
        save_path: str | None = None,
    ) -> None:
        """
        Print and/or save all matrices and stability results.

        Args:
            eigenvalues    : Output of compute_stability().
            is_stable      : Output of compute_stability().
            max_real       : Output of compute_stability().
            print_matrices : If True, print S, Phi, E, A, B to stdout.
            save_path      : If given, write all matrices to this txt file.
        """
        def _arr2str(a: np.ndarray) -> str:
            return np.array2string(
                np.asarray(a),
                max_line_width=99999,
                precision=6,
                suppress_small=True,
                threshold=99999,
            )

        def _section(title: str, content: str) -> str:
            bar = "=" * (len(title) + 4)
            return f"{bar}\n  {title}\n{bar}\n{content}\n"

        if print_matrices:
            for label, mat in [
                ("S",   self.S_H),
                ("Phi", self.Phi_H),
                ("E",   self.E),
                ("A",   self.A),
                ("B",   self.B),
            ]:
                if mat is not None:
                    print(f"\n--- {label} Matrix ---")
                    print(sparse.csc_matrix(mat))

            print("\n--- Stability Analysis ---")
            print(f"Is stable:     {is_stable}")
            print(f"Max Real Part: {max_real}")

        if save_path is not None:
            import os
            os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
            with open(save_path, "w", encoding="utf-8") as f:
                for label, mat in [
                    ("S Matrix",   self.S_H),
                    ("Phi Matrix", self.Phi_H),
                    ("E Matrix",   self.E),
                    ("A Matrix",   self.A),
                    ("B Matrix",   self.B),
                ]:
                    if mat is not None:
                        f.write(_section(label, _arr2str(mat)) + "\n")

                stability_txt = (
                    f"Is stable:     {is_stable}\n"
                    f"Max Real Part: {max_real}\n"
                )
                if eigenvalues is not None:
                    evs_txt = "\n".join(f"  {ev}" for ev in eigenvalues)
                    stability_txt += f"Eigenvalues:\n{evs_txt}\n"

                f.write(_section("Stability Analysis", stability_txt))

            print(f"\nMatrices saved to: {save_path}")

