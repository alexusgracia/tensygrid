"""
Unit tests for MatrixBuilder class.
Run with: pytest tests/test_matrix_builder.py -v
"""

import pytest
import sympy as sp
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.MatrixBuilder import MatrixBuilder


class TestMatrixBuilderBasic:
    """Tests for basic functionality of MatrixBuilder."""

    def test_simple_linear_equations(self):
        """Test with simple linear equations."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # Check symbols are correctly identified
        symbol_names = [str(s) for s in builder.symbols_list]
        assert 'dx1' in symbol_names
        assert 'dx2' in symbol_names
        assert 'x1' in symbol_names
        assert 'x2' in symbol_names
        assert 'u' in symbol_names
        
        # Check terms are correctly identified
        term_names = [str(t) for t in builder.terms_list]
        assert 'x1*x2' in term_names

    def test_matrix_S_dimensions(self):
        """Test that S matrix has correct dimensions."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        n_symbols = len(builder.symbols_list)
        n_terms = len(builder.terms_list)
        
        assert builder.S.shape == (n_symbols, n_terms)

    def test_matrix_P_dimensions(self):
        """Test that P matrix has correct dimensions."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        n_equations = len(equations)
        n_terms = len(builder.terms_list)
        
        assert builder.P.shape == (n_equations, n_terms)


class TestMatrixS:
    """Tests for S matrix construction."""

    def test_S_diagonal_ones(self):
        """Test that S matrix has 1s where symbols appear in their own term."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # For each basic symbol, it should have 1 in its own column
        for sym in builder.symbols_list:
            if sym in builder.term_to_index:
                i = builder.symbol_to_index[sym]
                j = builder.term_to_index[sym]
                assert builder.S[i, j] == 1

    def test_S_nonlinear_term_marking(self):
        """Test that S matrix correctly marks symbols in nonlinear terms."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # Find x1*x2 term index
        x1x2_term = None
        for term in builder.terms_list:
            if str(term) == 'x1*x2':
                x1x2_term = term
                break
        
        assert x1x2_term is not None
        j = builder.term_to_index[x1x2_term]
        
        # x1 should have 1 in x1*x2 column
        x1_sym = [s for s in builder.symbols_list if str(s) == 'x1'][0]
        i_x1 = builder.symbol_to_index[x1_sym]
        assert builder.S[i_x1, j] == 1
        
        # x2 should have 1 in x1*x2 column
        x2_sym = [s for s in builder.symbols_list if str(s) == 'x2'][0]
        i_x2 = builder.symbol_to_index[x2_sym]
        assert builder.S[i_x2, j] == 1

    def test_S_expected_matrix(self):
        """Test S matrix matches expected output."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # Expected S matrix (rows: dx1, dx2, x1, x2, u; cols: dx1, dx2, x1, x2, u, x1*x2)
        expected_S = sp.Matrix([
            [1, 0, 0, 0, 0, 0],  # dx1
            [0, 1, 0, 0, 0, 0],  # dx2
            [0, 0, 1, 0, 0, 1],  # x1
            [0, 0, 0, 1, 0, 1],  # x2
            [0, 0, 0, 0, 1, 0],  # u
        ])
        
        assert builder.S == expected_S


class TestMatrixP:
    """Tests for P matrix construction."""

    def test_P_coefficients_simple(self):
        """Test P matrix coefficients for simple equations."""
        equations = [
            '-dx1 - x1 + u',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # Expected P matrix
        expected_P = sp.Matrix([
            [-1, 0, -1, 0, 1, 0],   # -dx1 - x1 + u
            [0, -1, 0, -1, 0, 1],   # -dx2 - x2 + x1*x2
        ])
        
        assert builder.P == expected_P

    def test_P_coefficients_with_numeric(self):
        """Test P matrix with numeric coefficients."""
        equations = [
            '-dx1 - x1 + u + 2x2',
            '-dx2 - x2 + x1*x2',
        ]
        builder = MatrixBuilder(equations)
        
        # Expected P matrix
        expected_P = sp.Matrix([
            [-1, 0, -1, 2, 1, 0],   # -dx1 - x1 + u + 2*x2
            [0, -1, 0, -1, 0, 1],   # -dx2 - x2 + x1*x2
        ])
        
        assert builder.P == expected_P

    def test_P_negative_coefficients(self):
        """Test P matrix handles negative coefficients correctly."""
        equations = [
            '-3dx1 + 2x1 - 5u',
        ]
        builder = MatrixBuilder(equations)
        
        # Find column indices
        dx1_idx = builder.term_to_index[sp.Symbol('dx1')]
        x1_idx = builder.term_to_index[sp.Symbol('x1')]
        u_idx = builder.term_to_index[sp.Symbol('u')]
        
        assert builder.P[0, dx1_idx] == -3
        assert builder.P[0, x1_idx] == 2
        assert builder.P[0, u_idx] == -5


class TestImplicitMultiplication:
    """Tests for implicit multiplication preprocessing."""

    def test_implicit_multiplication_2x(self):
        """Test that 2x is converted to 2*x."""
        equations = ['2x1 + 3x2']
        builder = MatrixBuilder(equations)
        
        x1_idx = builder.term_to_index[sp.Symbol('x1')]
        x2_idx = builder.term_to_index[sp.Symbol('x2')]
        
        assert builder.P[0, x1_idx] == 2
        assert builder.P[0, x2_idx] == 3

    def test_implicit_multiplication_preserves_symbols(self):
        """Test that symbol names like dx1 are preserved."""
        equations = ['dx1 + x1']
        builder = MatrixBuilder(equations)
        
        symbol_names = [str(s) for s in builder.symbols_list]
        assert 'dx1' in symbol_names
        assert 'x1' in symbol_names


class TestSymbolOrdering:
    """Tests for symbol and term ordering."""

    def test_derivatives_first(self):
        """Test that derivative symbols (d...) come first."""
        equations = ['x1 + dx1 + u']
        builder = MatrixBuilder(equations)
        
        symbol_names = [str(s) for s in builder.symbols_list]
        dx1_idx = symbol_names.index('dx1')
        x1_idx = symbol_names.index('x1')
        u_idx = symbol_names.index('u')
        
        # Derivatives should come before states
        assert dx1_idx < x1_idx
        # States should come before inputs
        assert x1_idx < u_idx

    def test_nonlinear_terms_last(self):
        """Test that nonlinear terms come last in column ordering."""
        equations = ['x1*x2 + x1 + x2']
        builder = MatrixBuilder(equations)
        
        term_names = [str(t) for t in builder.terms_list]
        x1x2_idx = term_names.index('x1*x2')
        x1_idx = term_names.index('x1')
        x2_idx = term_names.index('x2')
        
        # Nonlinear term should come after basic symbols
        assert x1x2_idx > x1_idx
        assert x1x2_idx > x2_idx


class TestEquationWithEquals:
    """Tests for equations with '=' sign."""

    def test_equation_with_equals(self):
        """Test that equations with '=' are parsed correctly."""
        equations = ['dx1 = -x1 + u']
        builder = MatrixBuilder(equations)
        
        # This should be equivalent to 'dx1 + x1 - u = 0'
        dx1_idx = builder.term_to_index[sp.Symbol('dx1')]
        x1_idx = builder.term_to_index[sp.Symbol('x1')]
        u_idx = builder.term_to_index[sp.Symbol('u')]
        
        assert builder.P[0, dx1_idx] == 1
        assert builder.P[0, x1_idx] == 1
        assert builder.P[0, u_idx] == -1


class TestEdgeCases:
    """Tests for edge cases."""

    def test_empty_equations(self):
        """Test with empty equation list."""
        builder = MatrixBuilder([])
        assert builder.S.shape == (0, 0)
        assert builder.P.shape == (0, 0)

    def test_single_symbol(self):
        """Test with single symbol equation."""
        equations = ['x1']
        builder = MatrixBuilder(equations)
        
        assert len(builder.symbols_list) == 1
        assert len(builder.terms_list) == 1
        assert builder.P[0, 0] == 1

    def test_multiple_nonlinear_terms(self):
        """Test with multiple nonlinear terms."""
        equations = ['x1*x2 + x2*x3 + x1*x3']
        builder = MatrixBuilder(equations)
        
        # Should have 3 basic symbols
        assert len(builder.symbols_list) == 3
        # Only 3 terms (nonlinear only, basic symbols don't appear as separate terms)
        assert len(builder.terms_list) == 3
        
    def test_mixed_linear_nonlinear_terms(self):
        """Test with mixed linear and nonlinear terms."""
        equations = ['x1 + x2 + x1*x2']
        builder = MatrixBuilder(equations)
        
        # Should have 2 basic symbols
        assert len(builder.symbols_list) == 2
        # Should have 3 terms: x1, x2, x1*x2
        assert len(builder.terms_list) == 3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
