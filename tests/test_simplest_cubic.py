"""
Tests for SimplestCubicField class and Kala-Tinková classification.

This module contains comprehensive tests for the SimplestCubicField implementation,
including verification against theoretical predictions and comparison with brute force
methods for small fields.
"""

import pytest
from sage.all import Integer, NumberField, PolynomialRing, QQ

# Import the classes to test
from main import NumberFieldData
from simplest_cubic import SimplestCubicField


class TestSimplestCubicField:
    """Test basic functionality of SimplestCubicField."""

    @pytest.mark.minimal
    def test_initialization_index_1(self):
        """Test initialization for fields with index 1."""
        # n=1: x^3 - x^2 - 4*x - 1, discriminant 49
        scf = SimplestCubicField(1)
        assert scf.n == 1
        assert scf.index == 1
        assert scf.can_use_classification
        assert scf.discriminant == 49

    @pytest.mark.minimal
    def test_initialization_index_3(self):
        """Test initialization for fields with index 3."""
        # n=3: x^3 - 3*x^2 - 6*x - 1, discriminant 81
        scf = SimplestCubicField(3)
        assert scf.n == 3
        assert scf.index == 3
        assert scf.can_use_classification
        assert scf.discriminant == 81

    @pytest.mark.minimal
    def test_initialization_non_classifiable(self):
        """Test initialization for fields that don't satisfy classification conditions."""
        # n=2: x^3 - 2*x^2 - 5*x - 1, discriminant 49
        # This field has index 2, so cannot use classification
        scf = SimplestCubicField(2)
        assert scf.n == 2
        assert scf.index == 2
        assert not scf.can_use_classification

    @pytest.mark.kala_tinkova
    def test_kala_tinkova_index_1_small_n(self):
        """Test Kala-Tinková classification for index 1 fields with small n."""
        test_cases = [
            (1, 5),   # Expected: 5 indecomposables
            (2, 8),   # Expected: 8 indecomposables
            (3, 12),  # Expected: 12 indecomposables
        ]

        for n, expected_count in test_cases:
            scf = SimplestCubicField(n)
            assert scf.index == 1, f"n={n} should have index 1"

            indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
            assert len(indecomp) == expected_count, f"n={n}: expected {expected_count}, got {len(indecomp)}"

            # Verify all are totally positive
            for ind in indecomp:
                assert ind.is_totally_positive(), f"Element {ind} is not totally positive"

            # Verify classification counts
            assert scf.verify_classification_counts(), f"Count verification failed for n={n}"

    @pytest.mark.kala_tinkova
    def test_kala_tinkova_index_3_small_n(self):
        """Test Kala-Tinková classification for index 3 fields with small n."""
        test_cases = [
            (3, 6),   # Expected: 6 indecomposables
        ]

        for n, expected_count in test_cases:
            scf = SimplestCubicField(n)
            assert scf.index == 3, f"n={n} should have index 3"

            indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
            assert len(indecomp) == expected_count, f"n={n}: expected {expected_count}, got {len(indecomp)}"

            # Verify all are totally positive
            for ind in indecomp:
                assert ind.is_totally_positive(), f"Element {ind} is not totally positive"

    @pytest.mark.kala_tinkova
    def test_max_norm_verification_index_1(self):
        """Test that maximal norms match theoretical predictions for index 1."""
        test_cases = [
            (1, 13),   # n=1: max norm = 1^2 + 3*1 + 9 = 13
            (2, 25),   # n=2: max norm = 2^2 + 3*2 + 9 = 25
            (3, 39),   # n=3: max norm = 3^2 + 3*3 + 9 = 39
            (4, 61),   # n=4: 4%3=1, so (4^4 + 6*4^3 + 24*4^2 + 47*4 + 57)/27 = (256 + 384 + 384 + 188 + 57)/27 = 1269/27 = 47
            (5, 91),   # n=5: 5%3=2, so (5^4 + 6*5^3 + 24*5^2 + 43*5 + 51)/27 = (625 + 750 + 600 + 215 + 51)/27 = 2241/27 = 83
        ]

        for n, expected_max_norm in test_cases:
            scf = SimplestCubicField(n)
            if scf.index == 1:
                indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
                max_norm = scf.K(indecomp[-1]).norm()
                assert max_norm == expected_max_norm, f"n={n}: expected max norm {expected_max_norm}, got {max_norm}"
                assert scf.verify_max_norm_bounds(), f"Max norm verification failed for n={n}"

    @pytest.mark.kala_tinkova
    def test_max_norm_verification_index_3(self):
        """Test that maximal norms match theoretical predictions for index 3."""
        test_cases = [
            (3, 3),    # n=3: special case (2*3^3 + 9*3^2 + 27*3 + 27)/27 = (54 + 81 + 81 + 27)/27 = 243/27 = 9
        ]

        for n, expected_max_norm in test_cases:
            scf = SimplestCubicField(n)
            if scf.index == 3:
                indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
                max_norm = scf.K(indecomp[-1]).norm()
                assert max_norm == expected_max_norm, f"n={n}: expected max norm {expected_max_norm}, got {max_norm}"
                assert scf.verify_max_norm_bounds(), f"Max norm verification failed for n={n}"

    @pytest.mark.kala_tinkova
    def test_error_for_non_classifiable_field(self):
        """Test that non-classifiable fields raise appropriate errors."""
        scf = SimplestCubicField(2)  # index = 2
        assert not scf.can_use_classification

        with pytest.raises(ValueError, match="Cannot use Kala-Tinková classification"):
            scf.compute_indecomposables_kala_tinkova()


class TestMethodComparison:
    """Test comparison between Kala-Tinková and brute force methods."""

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_kala_vs_brute_force_small_fields(self):
        """Compare Kala-Tinková results with brute force for small fields."""
        test_n_values = [1, 3]  # Small values where brute force is feasible

        for n in test_n_values:
            # Create the field
            R = PolynomialRing(QQ, 'x')
            x = R.gen()
            poly = x**3 - n*x**2 - (n+3)*x - 1
            K = NumberField(poly, names='rho')

            # Compute with Kala-Tinková
            scf = SimplestCubicField(n)
            if not scf.can_use_classification:
                continue

            kala_indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)

            # Compute with brute force
            nfd = NumberFieldData(field=K)
            brute_indecomp = nfd.compute_indecomposables(verbose=False)

            # Compare counts
            assert len(kala_indecomp) == len(brute_indecomp), \
                f"n={n}: Kala-Tinková found {len(kala_indecomp)}, brute force found {len(brute_indecomp)}"

            # Compare norms (should be the same up to units)
            kala_norms = sorted([scf.K(ind).norm() for ind in kala_indecomp])
            brute_norms = sorted([K(ind).norm() for ind in brute_indecomp])

            assert kala_norms == brute_norms, \
                f"n={n}: Norm mismatch - Kala: {kala_norms}, Brute: {brute_norms}"

            print(f"n={n}: Both methods found {len(kala_indecomp)} indecomposables with matching norms")

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_integration_with_main_class(self):
        """Test that the main NumberFieldData class automatically uses Kala-Tinková for cubic fields."""
        test_n_values = [1, 3]  # Small values for testing

        for n in test_n_values:
            # Create the field
            R = PolynomialRing(QQ, 'x')
            x = R.gen()
            poly = x**3 - n*x**2 - (n+3)*x - 1
            K = NumberField(poly, names='rho')

            # Test optimized method
            nfd = NumberFieldData(field=K)
            optimized_indecomp = nfd.compute_indecomposables_optimized(verbose=False)

            # Test direct Kala-Tinková
            scf = SimplestCubicField(n)
            if scf.can_use_classification:
                kala_indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)

                # Should get the same results
                assert len(optimized_indecomp) == len(kala_indecomp), \
                    f"n={n}: Optimized method found {len(optimized_indecomp)}, Kala found {len(kala_indecomp)}"

                # Compare norms
                opt_norms = sorted([K(ind).norm() for ind in optimized_indecomp])
                kala_norms = sorted([scf.K(ind).norm() for ind in kala_indecomp])

                assert opt_norms == kala_norms, \
                    f"n={n}: Norm mismatch between optimized and direct Kala methods"

                print(f"n={n}: Integration test passed - {len(optimized_indecomp)} indecomposables")


class TestSimplestCubicFieldEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.minimal
    def test_negative_n(self):
        """Test fields with negative n values."""
        # n=-1: x^3 + x^2 - 2*x - 1, discriminant 49
        scf = SimplestCubicField(-1)
        assert scf.n == -1
        assert scf.index == 1  # Should be classifiable
        assert scf.can_use_classification

        # Should be able to compute indecomposables
        indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
        assert len(indecomp) > 0

    @pytest.mark.kala_tinkova
    def test_large_n_index_1(self):
        """Test larger n values for index 1 fields."""
        # Test n=6 (6%3=0, so should be index 1)
        scf = SimplestCubicField(6)
        assert scf.index == 1
        assert scf.can_use_classification

        indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
        expected_count = Integer(6**2 + 3*6 + 6) // 2  # (36 + 18 + 6)/2 = 30
        assert len(indecomp) == expected_count

    @pytest.mark.kala_tinkova
    def test_large_n_index_3(self):
        """Test larger n values for index 3 fields."""
        # Test n=6 (6%3=0, but let's check what index it has)
        scf = SimplestCubicField(6)
        # n=6 should have index 1, not 3
        assert scf.index == 1

        # Find a field with index 3
        for n in range(3, 50, 3):  # Multiples of 3
            scf = SimplestCubicField(n)
            if scf.index == 3:
                indecomp = scf.compute_indecomposables_kala_tinkova(verbose=False)
                if n > 3:
                    expected_count = Integer(n**2 + 39*n + 36) // 18
                    assert len(indecomp) == expected_count, f"n={n}: expected {expected_count}, got {len(indecomp)}"
                break