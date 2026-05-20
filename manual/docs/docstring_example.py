"""
docstring_example.py — Exemple de com documentar codi Python per a TenSyGrid
==============================================================================

Aquest fitxer serveix de referència per documentar correctament el codi.
La convenció del projecte és **NumPy-style docstrings**. Sphinx-AutoAPI
parseja aquests docstrings automàticament per generar la referència de l'API.

Estructura requerida:
 - Una línia de resum (una sola frase, verb en infinitiu, sense punt final)
 - Secció de descripció extesa (opcional)
 - Seccions estàndard: Parameters, Returns, Raises, Notes, Examples

Regles importants:
 - Usar raw strings (r\"\"\"...\"\"\") si el docstring conté LaTeX com \\frac, \\sum, etc.
 - Indentar el contingut del docstring al mateix nivell que les cometes.
 - No deixar línies amb indentació diferent a la base (trenca prepare_docstring).
"""

import numpy as np
from typing import Optional


# ---------------------------------------------------------------------------
# Exemple de funció senzilla
# ---------------------------------------------------------------------------

def compute_eigenvalues(A: np.ndarray, k: int = 10) -> tuple[np.ndarray, np.ndarray]:
    r"""Compute the k dominant eigenvalues of a matrix.

    Uses the implicitly restarted Arnoldi method (ARPACK) to find the ``k``
    eigenvalues with largest magnitude. Suitable for large sparse matrices
    where a full decomposition is impractical.

    Parameters
    ----------
    A : np.ndarray
        Square input matrix of shape ``(n, n)``. Can be dense or sparse.
    k : int, optional
        Number of eigenvalues to compute. Must satisfy ``k < n``.
        Default is ``10``.

    Returns
    -------
    eigenvalues : np.ndarray
        Complex array of shape ``(k,)`` with the dominant eigenvalues,
        sorted by descending magnitude.
    eigenvectors : np.ndarray
        Matrix of shape ``(n, k)`` whose columns are the right eigenvectors
        corresponding to each eigenvalue.

    Raises
    ------
    ValueError
        If ``k >= A.shape[0]`` or if ``A`` is not square.
    np.linalg.LinAlgError
        If ARPACK fails to converge within the maximum number of iterations.

    Notes
    -----
    The participation factor :math:`P_{ji}` of state :math:`i` in mode
    :math:`j` is defined as:

    .. math::

        P_{ji} = \\frac{w_{ji} \, v_{ji}}{w_j^H v_j}

    where :math:`v_j` and :math:`w_j` are the right and left eigenvectors
    of mode :math:`j`, respectively.

    Examples
    --------
    >>> import numpy as np
    >>> A = np.diag([-1.0, -2.0, -3.0, 0.5])
    >>> evals, evecs = compute_eigenvalues(A, k=2)
    >>> evals.shape
    (2,)
    """
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square, got shape {A.shape}")
    if k >= A.shape[0]:
        raise ValueError(f"k={k} must be less than matrix size {A.shape[0]}")
    eigenvalues, eigenvectors = np.linalg.eig(A)
    idx = np.argsort(np.abs(eigenvalues))[::-1][:k]
    return eigenvalues[idx], eigenvectors[:, idx]


# ---------------------------------------------------------------------------
# Exemple de classe
# ---------------------------------------------------------------------------

class StabilityAnalysis:
    r"""Perform small-signal stability analysis on a linearised power system.

    Wraps the eigenvalue computation and provides methods to assess system
    stability, compute participation factors, and export results.

    Parameters
    ----------
    A : np.ndarray
        State matrix of shape ``(n, n)``.
    E : np.ndarray, optional
        Descriptor matrix of shape ``(n, n)`` for generalised problems
        of the form :math:`E \dot{x} = A x`. If ``None``, the standard
        problem :math:`\dot{x} = A x` is solved.
    sigma : float, optional
        Shift-and-Invert spectral transformation centre. Eigenvalues near
        ``sigma`` are computed first. Default is ``0.0``.

    Attributes
    ----------
    n : int
        State-space dimension (number of state variables).
    is_stable : bool or None
        ``True`` if all computed eigenvalues have negative real part,
        ``False`` otherwise. ``None`` before :meth:`compute` is called.
    eigenvalues : np.ndarray or None
        Computed eigenvalues after calling :meth:`compute`.

    Examples
    --------
    >>> import numpy as np
    >>> A = np.diag([-1.0, -0.5, -3.0])
    >>> sa = StabilityAnalysis(A)
    >>> result = sa.compute(k=2)
    >>> sa.is_stable
    True
    """

    def __init__(
        self,
        A: np.ndarray,
        E: Optional[np.ndarray] = None,
        sigma: float = 0.0,
    ) -> None:
        self.A = A
        self.E = E
        self.sigma = sigma
        self.n = A.shape[0]
        self.is_stable: Optional[bool] = None
        self.eigenvalues: Optional[np.ndarray] = None

    def compute(self, k: int = 10) -> dict:
        r"""Compute eigenvalues and assess small-signal stability.

        Parameters
        ----------
        k : int, optional
            Number of eigenvalues to compute. Default is ``10``.

        Returns
        -------
        result : dict
            Dictionary with the following keys:

            ``eigenvalues`` : np.ndarray
                Complex array of shape ``(k,)``.
            ``stable`` : bool
                ``True`` if all real parts are strictly negative.
            ``margin`` : float
                Maximum real part (negative means stable).

        Notes
        -----
        Stability is assessed by the sign of the maximum real part:

        .. math::

            \\text{margin} = \\max_i \\, \\text{Re}(\\lambda_i)

        A system is small-signal stable if :math:`\\text{margin} < 0`.
        """
        evals, _ = compute_eigenvalues(self.A, k=k)
        self.eigenvalues = evals
        margin = float(np.max(evals.real))
        self.is_stable = margin < 0.0
        return {"eigenvalues": evals, "stable": self.is_stable, "margin": margin}

    def participation_factors(self) -> np.ndarray:
        r"""Compute the participation factor matrix.

        Participation factors quantify the contribution of each state variable
        to each oscillatory mode. For mode :math:`j` and state :math:`i`:

        .. math::

            P_{ji} = \\frac{w_{ji} \, v_{ji}}{w_j^H v_j}

        Returns
        -------
        P : np.ndarray
            Real matrix of shape ``(k, n)`` where ``k`` is the number of
            computed modes and ``n`` the number of state variables.

        Raises
        ------
        RuntimeError
            If :meth:`compute` has not been called yet.
        """
        if self.eigenvalues is None:
            raise RuntimeError("Call compute() before participation_factors().")
        evals, evecs_r = compute_eigenvalues(self.A, k=len(self.eigenvalues))
        _, evecs_l = compute_eigenvalues(self.A.T, k=len(self.eigenvalues))
        P = np.abs(evecs_l.conj() * evecs_r).T
        return P.real
