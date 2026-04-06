"""
Tests for the NumberFieldData class and main.py module.

These tests verify:
1. Number field construction and properties
2. Unit group computations (fundamental, totally positive)
3. Indecomposable enumeration
4. Data I/O formatting
"""

import pytest

import sys
import os

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

# Import Sage modules 
from sage.all import QQ, PolynomialRing, NumberField, ZZ

from main import NumberFieldData, create_field_from_polynomial, process_degree_file


class TestNumberFieldDataConstruction:
    """Test NumberFieldData initialization and basic properties."""
    
    #def test_init_with_label_only(self):
    #    """Test creating NumberFieldData with just an LMFDB label."""
    #    nfd = NumberFieldData(lmfdb_label="3.3.49.1")
    #    assert nfd.lmfdb_label == "3.3.49.1"
    #    assert nfd.K is None
    #    assert nfd.metadata == {}
    
    #def test_init_with_metadata(self):
    #    """Test creating NumberFieldData with metadata."""
    #    metadata = {
    #        'discriminant': 49,
    #        'regulator': 0.525,
    #        'class_number': 1,
    #        'degree': 3
    #    }
    #    nfd = NumberFieldData(lmfdb_label="3.3.49.1", metadata=metadata)
    #    assert nfd.discriminant == 49
    #    assert nfd.regulator == 0.525
    #    assert nfd.class_number == 1
    #    assert nfd.degree == 3
    
    def test_init_with_field(self):
        """Test creating NumberFieldData given coefficients"""
        # Real quadratic field with disc 5
        nfd = NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
        assert nfd.lmfdb_label == "2.2.5.1"
        assert nfd.K is not None
        assert nfd.degree == 2
        assert nfd.discriminant == 5
        assert nfd.class_number == 1
        assert abs(nfd.regulator - 0.48121182506) < 1e-6


class TestNumberFieldProperties:
    """Test number field properties."""
    
    @pytest.fixture
    def real_quadratic_field(self):
        """Fixture: Real quadratic field Q(sqrt(5))."""
        return NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
    
    @pytest.fixture
    def cubic_field(self):
        """Fixture: Totally real cubic field Q(a) where a^3 - a^2 - 2a + 1 = 0."""
        return NumberFieldData([1, -2, -1, 1], lmfdb_label="3.3.49.1")
    
    def test_discriminant_quadratic(self, real_quadratic_field):
        """Test discriminant computation for quadratic field."""
        disc = real_quadratic_field.discriminant
        assert disc == 5
    
    def test_discriminant_cubic(self, cubic_field):
        """Test discriminant computation for cubic field."""
        disc = cubic_field.discriminant
        assert disc == 49
    
    def test_degree_quadratic(self, real_quadratic_field):
        """Test degree for quadratic field."""
        assert real_quadratic_field.degree == 2
    
    def test_degree_cubic(self, cubic_field):
        """Test degree for cubic field."""
        assert cubic_field.degree == 3
    
    def test_regulator_quadratic(self, real_quadratic_field):
        """Test regulator computation for quadratic field."""
        reg = real_quadratic_field.regulator
        assert reg > 0
        # For Q(sqrt(5)), regulator should be log(phi) ~ 0.481
        assert 0.4 < reg < 0.6
    
    def test_class_number_quadratic(self, real_quadratic_field):
        """Test class number for quadratic field."""
        h = real_quadratic_field.class_number
        assert h == 1


class TestFundamentalUnits:
    """Test fundamental unit computations."""
    
    @pytest.fixture
    def real_quadratic_field(self):
        """Fixture: Real quadratic field Q(sqrt(5))."""
        return NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
    
    @pytest.fixture
    def cubic_field(self):
        """Fixture: Totally real cubic field Q(a) where a^3 - a^2 - 2a + 1 = 0."""
        return NumberFieldData([1, -2, -1, 1], lmfdb_label="3.3.49.1")
    
    def test_fundamental_units_quadratic(self, real_quadratic_field):
        """Test fundamental units for degree 2 field."""
        nfd = real_quadratic_field
        fun_units = nfd.fundamental_units
        assert len(fun_units) == 1  # degree - 1 = 2 - 1 = 1
        u = fun_units[0]
        # Verify it's a unit
        assert nfd.K(u).norm() in [1, -1]
    
    def test_fundamental_units_cubic(self, cubic_field):
        """Test fundamental units for degree 3 field."""
        nfd = cubic_field
        fun_units = nfd.fundamental_units
        assert len(fun_units) == 2  # degree - 1 = 3 - 1 = 2
        for u in fun_units:
            assert nfd.K(u).norm() in [1, -1]
    
    def test_fundamental_units_normalized_positive(self, real_quadratic_field):
        """Test that fundamental units are normalized to be positive."""
        nfd = real_quadratic_field
        fun_units = nfd.fundamental_units
        for u in fun_units:
            embeddings = nfd.K(u).complex_embeddings(nfd.precision)
            assert embeddings[0] > 0, f"Unit {u} not positive under first embedding"


class TestTotallyPositiveUnits:
    """Test totally positive unit computations."""
    
    @pytest.fixture
    def real_quadratic_field(self):
        """Fixture: Real quadratic field Q(sqrt(5))."""
        return NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
    
    @pytest.fixture
    def cubic_field(self):
        """Fixture: Totally real cubic field Q(a) where a^3 - a^2 - 2a + 1 = 0."""
        return NumberFieldData([1, -2, -1, 1], lmfdb_label="3.3.49.1")
    

    def test_tp_units_exist(self, real_quadratic_field):
        """Test that totally positive units are computed."""
        nfd = real_quadratic_field
        utp = nfd.totally_positive_units
        assert utp is not None
        assert len(utp) == nfd.degree - 1
    
    def test_tp_units_are_units(self, real_quadratic_field):
        """Test that totally positive units have norm 1."""
        nfd = real_quadratic_field
        utp = nfd.totally_positive_units
        for u in utp:
            assert nfd.K(u).norm() == 1
    
    def test_tp_units_are_totally_positive(self, real_quadratic_field):
        """Test that totally positive units are indeed totally positive."""
        nfd = real_quadratic_field
        utp = nfd.totally_positive_units
        for u in utp:
            embeddings = nfd.K(u).complex_embeddings(nfd.precision)
            assert all(e > 0 for e in embeddings)
            assert u.is_totally_positive()
    
    def test_tp_units_cubic(self, cubic_field):
        """Test totally positive units for cubic field."""
        nfd = cubic_field
        utp = nfd.totally_positive_units
        assert len(utp) == 2
        for u in utp:
            assert nfd.K(u).norm() == 1
            assert u.is_totally_positive()


class TestUnitRepresentatives:
    """Test unit representative computations."""
    
    @pytest.fixture
    def real_quadratic_field(self):
        """Fixture: Real quadratic field Q(sqrt(5))."""
        return NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
    
    def test_unit_representatives_exist(self, real_quadratic_field):
        """Test that unit representatives are computed."""
        nfd = real_quadratic_field
        reps = nfd.unit_representatives
        assert reps is not None
        assert len(reps) == 2**(nfd.degree)
    
    def test_unit_representatives_are_units(self, real_quadratic_field):
        """Test that all representatives are units."""
        nfd = real_quadratic_field
        reps = nfd.unit_representatives
        for rep in reps:
            norm = nfd.K(rep).norm()
            assert norm in [1, -1]
    
    def test_unit_representatives_distinct(self, real_quadratic_field):
        """Test that all representatives are distinct."""
        nfd = real_quadratic_field
        reps = nfd.unit_representatives
        assert len(reps) == len(set(str(r) for r in reps))


class TestDataFormatting:
    """Test data formatting and I/O."""
    
    @pytest.fixture
    def quadratic_with_indecomp(self):
        """Fixture: Quadratic field with computed indecomposables."""
        nfd = NumberFieldData([-5, 0, 1], lmfdb_label="2.2.5.1")
        # For testing, just set a simple indecomposable list
        nfd._indecomposables = [K(1), K(2)]
        return nfd
    
    def test_to_data_row_format(self, quadratic_with_indecomp):
        """Test formatting to data row."""
        nfd = quadratic_with_indecomp
        row = nfd.to_data_row()
        
        # Should be a string with pipe separators
        assert isinstance(row, str)
        assert "|" in row
        assert "2.2.5.1" in row
    
    def test_to_data_row_contains_required_fields(self, quadratic_with_indecomp):
        """Test that to_data_row contains key fields."""
        nfd = quadratic_with_indecomp
        row = nfd.to_data_row()
        
        # Should contain label, discriminant, regulator, class_number, indecomposables
        assert nfd.label in row
        assert "[" in row and "]" in row  # Indecomposables list


class TestCreateFieldFromPolynomial:
    """Test create_field_from_polynomial utility function."""
    
    def test_create_from_string_polynomial(self):
        """Test creating field from polynomial string."""
        nfd = create_field_from_polynomial('x**3 - 2')
        assert nfd.K is not None
        assert nfd.K.degree() == 3
    
    def test_create_from_sage_polynomial(self):
        """Test creating field from Sage polynomial object."""
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        poly = x**2 - 5
        
        nfd = create_field_from_polynomial(poly)
        assert nfd.K is not None
        assert nfd.K.degree() == 2


class TestIntegration:
    """Integration tests combining multiple components."""
    
    def test_full_workflow_quadratic(self):
        """Test complete workflow for quadratic field."""
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        K = QQ.extension(x**2 - 5, names='a')
        
        nfd = NumberFieldData(label="2.2.5.1", field=K)
        
        # Check all properties can be accessed
        assert nfd.label == "2.2.5.1"
        assert nfd.degree == 2
        assert nfd.discriminant == 5
        assert nfd.regulator > 0
        assert nfd.class_number >= 1
        assert nfd.fundamental_units is not None
        assert nfd.totally_positive_units is not None
        assert nfd.unit_representatives is not None
    
    def test_configuration_options(self):
        """Test that configuration options can be set."""
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        K = QQ.extension(x**2 - 5, names='a')
        
        nfd = NumberFieldData(field=K)
        nfd.exit_for_nonunit = True
        nfd.check_asserts = True
        nfd.precision = 200
        
        assert nfd.exit_for_nonunit == True
        assert nfd.check_asserts == True
        assert nfd.precision == 200


# ============================================================================
# Minimal/Fast Tests (for CI)
# ============================================================================

class TestMinimal:
    """Minimal tests that run quickly on CI without heavy Sage computations."""
    
    def test_number_field_data_instantiation(self):
        """Test that NumberFieldData can be instantiated."""
        nfd = NumberFieldData(label="test")
        assert nfd.label == "test"
    
    def test_import_main_module(self):
        """Test that main module imports without errors."""
        import main
        assert hasattr(main, 'NumberFieldData')
        assert hasattr(main, 'create_field_from_polynomial')
        assert hasattr(main, 'batch_compute_indecomposables')
    
    def test_sage_availability(self):
        """Test Sage availability and graceful degradation."""
        try:
            from sage.all import QQ, ZZ, Integer
            SAGE_AVAILABLE = True
            
            # Test basic Sage functionality if available
            assert QQ(1/2) == QQ(1)/QQ(2)
            assert ZZ(5).is_prime()
            assert Integer(10) == 10
            
        except ImportError:
            SAGE_AVAILABLE = False
        
        # This test should pass whether Sage is available or not
        # It just checks that the import logic works
        assert isinstance(SAGE_AVAILABLE, bool)
 

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
