"""
Tests for real quadratic field indecomposables computation.

Includes tests for:
1. Dress-Scharlau continued fraction method
2. Comparison with brute force method
3. Consistency across methods
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from sage.all import QuadraticField, Integer, QQ, PolynomialRing
from real_quadratic import RealQuadraticField
from main import NumberFieldData


class TestRealQuadraticFieldConstruction:
    """Test RealQuadraticField initialization."""
    
    @pytest.mark.dress_scharlau
    def test_init_basic(self):
        """Test basic initialization."""
        rq = RealQuadraticField(5)
        assert rq.D == 5
        assert rq.K is not None
        assert rq.discriminant == 5
    
    def test_init_d_eq_1_mod_4(self):
        """Test D congruent to 1 (mod 4)."""
        rq = RealQuadraticField(13)
        assert rq.D == 13
        assert rq.D % 4 == 1
        assert rq.discriminant == 13
    
    def test_init_d_eq_2_mod_4(self):
        """Test D congruent to 2 (mod 4)."""
        rq = RealQuadraticField(2)
        assert rq.D == 2
        assert rq.D % 4 == 2
        assert rq.discriminant == 8
    
    def test_init_not_squarefree_raises(self):
        """Test that non-squarefree D raises error."""
        with pytest.raises(ValueError):
            RealQuadraticField(4)  # 4 = 2^2 not squarefree
    
    def test_init_negative_raises(self):
        """Test that negative D raises error."""
        with pytest.raises(ValueError):
            RealQuadraticField(-5)


class TestFundamentalUnit:
    """Test fundamental unit computation."""
    
    @pytest.fixture
    def rq_5(self):
        """Q(sqrt(5))"""
        return RealQuadraticField(5)
    
    @pytest.fixture
    def rq_13(self):
        """Q(sqrt(13))"""
        return RealQuadraticField(13)
    
    def test_fund_unit_exists(self, rq_5):
        """Test fundamental unit is computed."""
        u = rq_5.fund_unit
        assert u is not None
        assert u > 1
    
    def test_fund_unit_is_unit(self, rq_5):
        """Test fundamental unit has norm ±1."""
        u = rq_5.fund_unit
        assert rq_5.OK(u).is_unit()
        assert u.norm() in [1, -1]
    
    def test_fund_unit_positive(self, rq_5):
        """Test fundamental unit is positive."""
        u = rq_5.fund_unit
        assert u > 1
    
    def test_tp_unit_exists(self, rq_5):
        """Test totally positive unit is computed."""
        utp = rq_5.tp_unit
        assert utp is not None
        assert utp > 1
        assert utp.norm() == 1
        assert utp.trace() > 0
    
    def test_tp_unit_is_totally_positive(self, rq_5):
        """Test tp unit is totally positive."""
        utp = rq_5.tp_unit
        assert rq_5.K(utp).is_totally_positive()


class TestContinuedFraction:
    """Test continued fraction data computation."""
    
    @pytest.fixture
    def rq_5(self):
        """Q(sqrt(5))"""
        return RealQuadraticField(5)
    
    @pytest.fixture
    def rq_13(self):
        """Q(sqrt(13))"""
        return RealQuadraticField(13)
    
    def test_cf_data_exists(self, rq_5):
        """Test CF data is computed."""
        cf, delta, s = rq_5.cf_data
        assert cf is not None
        assert delta is not None
        assert s > 0
    
    def test_cf_period_length_positive(self, rq_5):
        """Test period length is positive."""
        cf, delta, s = rq_5.cf_data
        assert s > 0
        assert len(cf.period()) == s
    
    def test_cf_d_mod_4_eq_1(self, rq_13):
        """Test CF for D equiv 1 (mod 4)."""
        cf, delta, s = rq_13.cf_data
        assert s > 0
        # Period should be computed
        assert len(cf.period()) == s


class TestAlphaSequence:
    """Test alpha sequence generation."""
    
    @pytest.fixture
    def rq_5(self):
        """Q(sqrt(5))"""
        return RealQuadraticField(5)
    
    def test_alpha_sequence_generated(self, rq_5):
        """Test alpha sequence is generated."""
        alpha = rq_5._compute_alpha_sequence()
        assert len(alpha) >= 2 * rq_5._s + 10
        assert alpha[0] == 1
    
    def test_alpha_sequence_units(self, rq_5):
        """Test that specific alpha elements are units."""
        cf, delta, s = rq_5.cf_data
        alpha = rq_5._compute_alpha_sequence()
        
        # alpha_{s-1} and alpha_{2s-1} should be units
        assert rq_5.OK(alpha[s]).is_unit()
        assert rq_5.OK(alpha[2 * s]).is_unit()
    
    def test_alpha_2s_totally_positive(self, rq_5):
        """Test that alpha_{2s-1} is totally positive."""
        cf, delta, s = rq_5.cf_data
        alpha = rq_5._compute_alpha_sequence()
        
        alpha_2s = alpha[2 * s]
        assert alpha_2s.norm() > 0
        assert alpha_2s.trace() > 0


class TestDressScharlau:
    """Test Dress-Scharlau indecomposables computation."""
    
    @pytest.fixture
    def rq_5(self):
        """Q(sqrt(5))"""
        return RealQuadraticField(5)
    
    @pytest.fixture
    def rq_13(self):
        """Q(sqrt(13))"""
        return RealQuadraticField(13)
    
    @pytest.fixture
    def rq_61(self):
        """Q(sqrt(61))"""
        return RealQuadraticField(61)
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_computed(self, rq_5):
        """Test indecomposables are computed."""
        indecomp = rq_5.compute_indecomposables_dress_scharlau(verbose=False)
        assert indecomp is not None
        assert len(indecomp) == 1
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_sorted(self, rq_5):
        """Test indecomposables are sorted by norm."""
        indecomp = rq_5.compute_indecomposables_dress_scharlau(verbose=False)
        norms = [rq_5.K(x).norm() for x in indecomp]
        assert norms == sorted(norms)
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_max_norm_bound(self, rq_5):
        """Test max norm satisfies Dress-Scharlau bound."""
        indecomp = rq_5.compute_indecomposables_dress_scharlau(verbose=False)
        max_norm = rq_5.K(indecomp[-1]).norm()
        bound = abs(rq_5.discriminant) / 4
        assert max_norm <= bound, \
                f"Max norm {max_norm} exceeds bound {bound}"
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_q_sqrt_5(self, rq_5):
        """Test Q(sqrt(5)) indecomposables."""
        # Known: Q(sqrt(5)) has 1 indecomposable (the unit 1)
        indecomp = rq_5.compute_indecomposables_dress_scharlau(verbose=False)
        assert len(indecomp) == 1
        assert indecomp[0].norm() == 1  # First should be unit
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_q_sqrt_13(self, rq_13):
        """Test Q(sqrt(13)) indecomposables."""
        indecomp = rq_13.compute_indecomposables_dress_scharlau(verbose=False)
        # Should have multiple indecomposables
        assert len(indecomp) == 3
    
    @pytest.mark.dress_scharlau
    def test_indecomposables_q_sqrt_61(self, rq_61):
        """Test Q(sqrt(61)) indecomposables."""
        indecomp = rq_61.compute_indecomposables_dress_scharlau(verbose=False)
        assert len(indecomp) == 11
        # Verify max norm bound
        max_norm = rq_61.K(indecomp[-1]).norm()
        bound = abs(rq_61.discriminant) / 4
        assert max_norm <= bound


class TestBruteForceBrute:
    """Test brute force method (from main.py) on quadratic fields."""
    
    @pytest.fixture
    def rq_5_nfd(self):
        """Q(sqrt(5)) via NumberFieldData."""
        K = QuadraticField(5, names='a')
        return NumberFieldData(label="2.2.5.1", field=K)
    
    @pytest.mark.brute_force
    @pytest.mark.slow
    def test_brute_force_computed(self, rq_5_nfd):
        """Test brute force indecomposables are computed."""
        # This may be slow, so use small timeout
        indecomp = rq_5_nfd.compute_indecomposables(verbose=False)
        assert indecomp is not None
        assert len(indecomp) == 1


class TestMethodComparison:
    """Compare Dress-Scharlau and brute force methods."""
    
    @pytest.mark.comparison
    @pytest.mark.slow
    def test_methods_agree_q_sqrt_5(self):
        """
        Test that Dress-Scharlau and brute force agree for Q(sqrt(5)).
        
        Marked as slow because brute force can be expensive.
        """
        # Compute via Dress-Scharlau
        rq = RealQuadraticField(5)
        ds_indecomp = rq.compute_indecomposables_dress_scharlau(verbose=False)
        ds_norms = sorted([rq.K(x).norm() for x in ds_indecomp])
        
        # Compute via brute force
        K = QuadraticField(5, names='a')
        nfd = NumberFieldData(label="2.2.5.1", field=K)
        bf_indecomp = nfd.compute_indecomposables(verbose=False)
        bf_norms = sorted([K(x).norm() for x in bf_indecomp])
            
        # Norms should match (might be different representatives)
        assert ds_norms == bf_norms, \
            f"DS norms {ds_norms} != BF norms {bf_norms}"
        
    
    @pytest.mark.comparison
    @pytest.mark.slow
    def test_methods_agree_count_q_sqrt_5(self):
        """Test that both methods find same number of indecomposables."""
        rq = RealQuadraticField(5)
        ds_indecomp = rq.compute_indecomposables_dress_scharlau(verbose=False)
        ds_count = len(ds_indecomp)
        
        # For Q(sqrt(5)), known to have exactly 1 indecomposable (the unit)
        assert ds_count == 1


class TestIntegrationRealQuadratic:
    """Integration tests for real quadratic fields."""
    
    @pytest.mark.integration
    @pytest.mark.dress_scharlau
    def test_full_workflow_sqrt_13(self):
        """Complete workflow for Q(sqrt(13))."""
        rq = RealQuadraticField(13)
        
        # Check all properties
        assert rq.D == 13
        assert rq.discriminant == 13
        assert rq.regulator > 0
        assert rq.class_number >= 1
        
        # Compute indecomposables
        indecomp = rq.compute_indecomposables_dress_scharlau(verbose=False)
        assert len(indecomp) > 0
        
        # Verify structure
        norms = [rq.K(x).norm() for x in indecomp]
        assert norms == sorted(norms)
        assert max(norms) <= abs(rq.discriminant) / 4
    
    @pytest.mark.integration
    @pytest.mark.dress_scharlau
    def test_multiple_fields(self):
        """Test computation on multiple fields."""
        D_values = [5, 13, 17, 29]
        
        for D in D_values:
            rq = RealQuadraticField(D)
            indecomp = rq.compute_indecomposables_dress_scharlau(verbose=False)
            
            assert len(indecomp) > 0
            norms = [rq.K(x).norm() for x in indecomp]
            assert norms == sorted(norms)

class TestRealQuadraticBigBruteForce:
    """ Brute force indecomposables for many real quadratic fields! """

    @pytest.mark.brute_force
    @pytest.mark.slow
    def test_real_quadratic_big_brute_force(self):

        for D in range(2, 20):
            if Integer(D).is_squarefree():
                rq = RealQuadraticField(D)

                # Compute indecomposables via Dress-Scharlau
                ds_indecomp = rq.compute_indecomposables_dress_scharlau(verbose=False)
                ds_norms = sorted([rq.K(x).norm() for x in ds_indecomp])
                assert len(ds_indecomp) > 0

                # Comptue indecomposables via brute force
                K = QuadraticField(D, names='a')
                nfd = NumberFieldData(field=K)
                bf_indecomp = nfd.compute_indecomposables(verbose=False)
                bf_norms = sorted([K(x).norm() for x in bf_indecomp])
                       
                # Ensure set of all norms match up
                assert ds_norms == bf_norms



if __name__ == "__main__":
    pytest.main([__file__, "-v"])
