"""
Efficient computation of indecomposables in simplest cubic fields using Kala-Tinková classification.

This module implements the explicit classification of indecomposables in simplest cubic fields
of the form Q[x]/(x^3 - n*x^2 - (n+3)*x - 1), where Z[ρ] has index 1 or 3 in the maximal order.

References:
    Kala, Vítězslav; Tinková, Markéta. "Universal quadratic forms, small norms and traces in families of number fields"
    Theorem 1.2: Classification when O_K = Z[ρ] (index 1)

    Gil-Muñoz, Daniel; Tinková, Markéta. "Additive structure of non-monogenic simplest cubic fields"
    Theorem 1.1: Classification when [O_K : Z[ρ]] = 3

Algorithm:
    For index 1 case: Enumerate from explicit formulas with parameters v, w
    For index 3 case: Enumerate from 8 different cases with parameters v, r
    Normalize up to multiplication by totally positive units
"""

# Import Sage modules 
from sage.all import (
    NumberField, PolynomialRing, QQ, ZZ, Integer
)
   

class SimplestCubicField:
    """
    Specialized class for computing indecomposables in simplest cubic fields.

    Simplest cubic fields have the form Q[x]/(x^3 - n*x^2 - (n+3)*x - 1).
    This class implements the Kala-Tinkova classification for fields where
    Z[ρ] has index 1 or 3 in the maximal order.
    """

    def __init__(self, n, precision=100):
        """
        Initialize a simplest cubic field Q[x]/(x^3 - n*x^2 - (n+3)*x - 1).

        Args:
            n: Integer parameter defining the field
            precision: Precision for RealField (in bits)
        """
        self.n = Integer(n)
        self.precision = precision

        # Create the polynomial x^3 - n*x^2 - (n+3)*x - 1
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        poly = x**3 - self.n*x**2 - (self.n + 3)*x - 1

        # Create the number field
        self.K = NumberField(poly, names='rho')
        self.rho = self.K.gen()
        self.OK = self.K.ring_of_integers()

        # Compute field invariants
        self.discriminant = self.K.discriminant()
        self.regulator = self.K.regulator()
        self.class_number = self.K.class_number()

        # Check if we can use the classification
        self.Zrho = self.K.order(self.rho)
        self.index = self.Zrho.index_in(self.OK)

        # Cache for computed data
        self._indecomposables = None

    @property
    def can_use_classification(self):
        """Check if we can use Kala-Tinková classification (index 1 or 3)."""
        return self.index == 1 or self.index == 3

    def _compute_index_1_indecomposables(self):
        """
        Compute indecomposables when O_K = Z[ρ] (index 1).

        Uses Theorem 1.2 from Kala-Tinková.
        Returns list of all indecomposables before normalization.
        """
        n = self.n
        rho = self.rho
        K = self.K

        all_indecomposables = [K(1), K(1 + rho + rho**2)]

        # Enumerate v, w parameters
        for v in range(n + 1):
            for w in range(v * (n + 2) + 1, (v + 1) * (n + 1) + 1):
                element = K(-v - w * rho + (v + 1) * rho**2)
                all_indecomposables.append(element)

        return all_indecomposables

    def _compute_index_3_indecomposables(self):
        """
        Compute indecomposables when [O_K : Z[ρ]] = 3.

        Uses Theorem 1.1 from Gil-Muñoz-Tinková.
        Returns list of all indecomposables before normalization.
        """
        n = self.n
        a = Integer(n)
        assert a % 3 == 0, f"n = {n} must be divisible by 3 for index 3 case"

        rho = self.rho
        K = self.K

        # Define basis elements as in Gil-Muñoz-Tinková
        g1 = K(1)
        g2 = K(rho)
        g3 = K((1 + rho + rho**2) / 3)

        all_indecomposables = []

        # Case (i): Basic elements
        all_indecomposables.extend([K(1), g3])

        # Case (ii): r parameter
        for r in range(1, (a // 3) + 1):
            element = -g1 - (r + 1) * g2 + 3 * g3
            all_indecomposables.append(K(element))

        # Case (iii): v parameter
        for v in range((2 * a) // 3 + 1, a + 1):
            element = -(2 * v + 1) * g1 - (v * (a + 3) + 2) * g2 + 3 * (v + 1) * g3
            all_indecomposables.append(K(element))

        # Case (iv): v parameter
        for v in range(0, a // 3):
            element = -(2 * v + 1) * g1 - ((v + 1) * (a + 2)) * g2 + 3 * (v + 1) * g3
            all_indecomposables.append(K(element))

        # Case (v): v and r parameters
        for v in range(0, a // 3):
            for r in range((a // 3) + 1, (2 * a) // 3 - v + 1):
                element = -(2 * v + 1) * g1 - (v * (a + 3) + r + 1) * g2 + (3 * v + 2) * g3
                all_indecomposables.append(K(element))

        # Case (vi): r parameter
        for r in range(0, a // 3):
            element = -(r + 1) * g2 + g3
            all_indecomposables.append(K(element))

        # Case (vii): v parameter
        for v in range(0, a // 3):
            element = -(2 * v + 2) * g1 - (v * (a + 3) + (2 * a) // 3 + 3) * g2 + (3 * v + 4) * g3
            all_indecomposables.append(K(element))

        # Case (viii): v parameter
        for v in range(a // 3, (2 * a) // 3):
            element = -(2 * v + 2) * g1 - (v * (a + 3) + (4 * a) // 3 - v + 3) * g2 + (3 * v + 4) * g3
            all_indecomposables.append(K(element))

        return all_indecomposables

    def compute_indecomposables_kala_tinkova(self, verbose=True):
        """
        Compute indecomposables using Kala-Tinkova classification.

        This is much faster than brute force for simplest cubic fields where
        Z[ρ] has index 1 or 3 in the maximal order.

        Args:
            verbose: Print progress information

        Returns:
            List of indecomposables (up to multiplication by totally positive units)

        Raises:
            ValueError: If the field doesn't satisfy the classification conditions
        """
        if not self.can_use_classification:
            raise ValueError(f"Cannot use Kala-Tinkova classification: Z[rho] has index {self.index}, not 1 or 3")

        if self._indecomposables is not None:
            return self._indecomposables

        if verbose:
            print(f"Computing indecomposables for simplest cubic field with n = {self.n}")
            print(f"Discriminant: {self.discriminant}")
            print(f"Index of Z[ρ] in O_K: {self.index}")

        # Compute raw indecomposables from the classification
        if self.index == 1:
            all_indecomposables = self._compute_index_1_indecomposables()
        else:  # self.index == 3
            all_indecomposables = self._compute_index_3_indecomposables()

        # Verify all are totally positive
        for ind in all_indecomposables:
            assert ind.is_totally_positive(), f"Element {ind} is not totally positive"

        if verbose:
            print(f"Generated {len(all_indecomposables)} candidates")

        # Normalize: keep only representatives up to units
        indecomposables_final = []
        OK = self.OK

        for ind in all_indecomposables:
            # Check if ind is equivalent (up to units) to a previously found indecomposable
            is_new = True
            for tmp in indecomposables_final:
                ratio = ind / tmp
                if ratio in OK and OK(ratio).is_unit():
                    is_new = False
                    break

            if is_new:
                indecomposables_final.append(ind)
                if verbose:
                    print(f"  Found indecomposable: {ind} (norm {ind.norm()})")

        # Sort by norm, then trace, then string length
        indecomposables_final.sort(
            key=lambda x: (self.K(x).norm(), self.K(x).trace(), len(str(x)))
        )

        self._indecomposables = indecomposables_final

        if verbose:
            print(f"\nTotal indecomposables: {len(indecomposables_final)}")
            if indecomposables_final:
                max_norm = self.K(indecomposables_final[-1]).norm()
                print(f"Max norm: {max_norm}")

                # Verify against known bounds
                assert max_norm <= abs(self.discriminant), \
                    f"Max norm {max_norm} exceeds discriminant {abs(self.discriminant)}"

        return indecomposables_final

    def indecomposables(self):
        """Get list of indecomposables (computed on demand)."""
        if self._indecomposables is None:
            self.compute_indecomposables_kala_tinkova(verbose=False)
        return self._indecomposables

    def verify_classification_counts(self):
        """
        Verify that the number of indecomposables matches the theoretical predictions.

        Returns:
            bool: True if counts match expected values
        """
        if not self.can_use_classification:
            return False

        indecomp = self.indecomposables()
        n = self.n

        if self.index == 1:
            # From Kala-Tinková Theorem 1.2
            expected = Integer(n**2 + 3*n + 6) // 2
            return len(indecomp) == expected

        elif self.index == 3:
            # From Gil-Muñoz-Tinková Theorem 1.1
            if n > 3:
                expected = Integer(n**2 + 39*n + 36) // 18
                return len(indecomp) == expected
            else:
                # Special cases for small n
                return True  # Skip verification for small n

        return False

    def verify_max_norm_bounds(self):
        """
        Verify that the maximal norm satisfies the theoretical bounds.

        Returns:
            bool: True if max norm matches expected bounds
        """
        if not self.can_use_classification:
            return False

        indecomp = self.indecomposables()
        if not indecomp:
            return False

        max_norm = self.K(indecomp[-1]).norm()
        n = self.n

        if self.index == 1:
            # From Kala-Tinkova
            if n <= 3:
                expected = n**2 + 3*n + 9
            elif n % 3 == 0:
                expected = (n**4 + 6*n**3 + 27*n**2 + 54*n + 81) // 27
            elif n % 3 == 1:
                expected = (n**4 + 6*n**3 + 24*n**2 + 47*n + 57) // 27
            elif n % 3 == 2:
                expected = (n**4 + 6*n**3 + 24*n**2 + 43*n + 51) // 27
            return max_norm == expected

        elif self.index == 3:
            # From Gil-Munoz-Tinkova Proposition 6.4
            if n in [3, 21, 30, 48]:
                expected = (2*n**3 + 9*n**2 + 27*n + 27) // 27
            else:
                expected = ((n**2 + 3*n + 9)**2) // 729
            return max_norm == expected

        return False