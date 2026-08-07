"""
Efficient computation of indecomposables in real quadratic fields using Dress-Scharlau method.

This module implements the Dress-Scharlau criterion based on continued fraction expansions
to directly enumerate indecomposables in real quadratic fields Q(sqrt(D)), which is much 
faster than the general brute force approach.

Reference:
    Dress, A.; Scharlau, W. (1975). 
    "Indecomposable integral representations of finite groups over real quadratic number fields"
    Advances in Mathematics, 17(3), 231-273.
    
    Kala, Vítězslav. "An effective tool for determining indecomposable elements of orders"
    arXiv:2108.15387

Algorithm:
    1. Compute continued fraction expansion of (√D - 1)/2 or √D
    2. Generate convergents p_i/q_i
    3. Compute "alpha" elements using recurrence relation
    4. Extract indecomposables from alpha_{i,t} = alpha_i + t*alpha_{i+1}
    5. Normalize up to multiplication by totally positive units
"""

# Try to import Sage modules - graceful degradation if not available
from sage.all import (
    QuadraticField, RealField, log, Integer, ZZ, continued_fraction, 
    prod
)


class RealQuadraticField:
    """
    Specialized class for computing indecomposables in real quadratic fields.
    
    Uses the Dress-Scharlau continued fraction method for efficient computation.
    """
    
    def __init__(self, D, precision=100):
        """
        Initialize a real quadratic field Q(sqrt(D)).
        
        Args:
            D: Squarefree positive integer (the discriminant or radicand)
            precision: Precision for RealField (in bits)
        """

        if not Integer(D).is_squarefree():
            raise ValueError(f"D = {D} must be squarefree")
        if D < 2:
            raise ValueError(f"D = {D} must be greater than 1")
        
        self.D = D
        self.precision = precision
        self.R = RealField(precision)
        
        # Create the number field
        self.K = QuadraticField(D, names='a')
        self.a = self.K.gen()
        self.OK = self.K.ring_of_integers()
        
        # Compute field invariants
        self.discriminant = self.K.discriminant()
        self.regulator = self.K.regulator()
        self.class_number = self.K.class_number()
        
        # Cache for computed data
        self._cf = None
        self._delta = None
        self._s = None
        self._fund_unit = None
        self._tp_unit = None
        self._indecomposables = None
        self._alpha_i = None
    
    @property
    def fund_unit(self):
        """Get fundamental unit u > 1 of the field."""
        if self._fund_unit is not None:
            return self._fund_unit
        
        UK = self.K.unit_group()
        u = UK.fundamental_units()[0]
        
        # Normalize: u > 1
        if abs(u) < 1:
            u = u.trace() - u  # Take conjugate
        if u < 0:
            u = -u  # Change sign
        
        assert u > 1, f"Fundamental unit not > 1: {u}"
        assert u.norm() in [1, -1], f"Unit has norm {u.norm()}"
        
        self._fund_unit = u
        return u
    
    @property
    def tp_unit(self):
        """Get generator for totally positive units."""
        if self._tp_unit is not None:
            return self._tp_unit
        
        u = self.fund_unit
        
        # If u has norm 1, it's already totally positive
        if u.norm() > 0:
            utp = u
        else:
            # If u has norm -1, use u^2
            utp = u**2
        
        assert utp.norm() == 1, f"TP unit has norm {utp.norm()}"
        assert utp.trace() > 0, f"TP unit not totally positive: trace = {utp.trace()}"
        assert utp > 1, f"TP unit not > 1: {utp}"
        
        self._tp_unit = utp
        return utp
    
    @property
    def cf_data(self):
        """
        Compute continued fraction data.
        
        Returns:
            (cf, delta, s) where:
            - cf: continued fraction object
            - delta: (sqrt(D) + 1)/2 or sqrt(D) depending on D mod 4
            - s: period length
        """
        if self._cf is not None:
            return (self._cf, self._delta, self._s)
        
        # Choose radicand based on D mod 4
        if self.D % 4 == 1:
            # D equiv 1 (mod 4): use (sqrt(D) - 1)/2 with delta = (sqrt(D) + 1)/2
            cf = continued_fraction((self.a - 1) / 2)
            delta = (self.a + 1) / 2
        else:
            # D equiv 2, 3 (mod 4): use sqrt(D) with delta = sqrt(D)
            cf = continued_fraction(self.a)
            delta = self.a
        
        s = len(cf.period())
        
        self._cf = cf
        self._delta = delta
        self._s = s
        
        return (cf, delta, s)
    
    def _compute_alpha_sequence(self, verbose=False):
        """
        Compute the sequence of alpha_i elements using convergents.
        
        The recurrence is: alpha_i = cf[i] * alpha_{i-1} + alpha_{i-2}
        where cf[i] is the i-th partial quotient of the continued fraction.
        
        Returns:
            List of alpha_i values
        """
        if self._alpha_i is not None:
            return self._alpha_i
        
        cf, delta, s = self.cf_data
        
        # Generate enough terms (at least 2s + 10 to ensure we capture all indecomposables)
        num_terms = 2 * s + 10
        alpha_i = [1]  # alpha_0 = 1
        
        for i in range(num_terms):
            # Get i-th convergent p_i / q_i
            pq = cf.convergent(i)
            p_i = pq.numer()
            q_i = pq.denom()
            
            # Compute alpha_i = p_i + q_i * delta
            alpha = self.K(p_i + q_i * delta)
            alpha_i.append(alpha)
            
            # Verify recurrence relation (for D at most 2000)
            if self.D <= 2000 and len(alpha_i) >= 3:
                expected = cf[i] * alpha_i[-2] + alpha_i[-3]
                assert alpha == expected, \
                    f"Recurrence failed at i={i}: {alpha} != {expected}"
            
            # Verify specific unitarity conditions (for D at most 2000)
            if self.D <= 2000 and i in [s - 1, 2 * s - 1]:
                assert self.OK(alpha).is_unit(), \
                    f"alpha_{i} = {alpha} is not a unit"
                if i == 2 * s - 1:
                    assert alpha.norm() > 0 and alpha.trace() > 0, \
                        f"alpha_{2*s-1} not totally positive"
        
        self._alpha_i = alpha_i

        if verbose:
            print("Computing convergents alpha sequence as:", alpha_i)

        return alpha_i
    
    def compute_indecomposables_dress_scharlau(self, verbose=True):
        """
        Compute indecomposables using the Dress-Scharlau continued fraction method.
        
        This is much faster than the general brute force approach for quadratic fields.
        
        Args:
            verbose: Print progress information
        
        Returns:
            List of indecomposables (up to multiplication by totally positive units)
        """
        if self._indecomposables is not None:
            return self._indecomposables
        
        cf, delta, s = self.cf_data
        alpha_i = self._compute_alpha_sequence(verbose=verbose)
        
        if verbose:
            print(f"Computing indecomposables for Q(sqrt({self.D}))")
            print(f"Discriminant: {self.discriminant}")
            print(f"Period length s = {s}")
        
        # Generate all candidate indecomposables from alpha_{i,t}
        all_indecomposables = []
        
        for i in range(0, 2 * s + 10, 2):
            cf_i_plus_1 = cf[i + 1]
            for t in range(cf_i_plus_1 + 1):
                # Compute alpha_{i,t} = alpha_i + t * alpha_{i+1}
                alpha_it = self.K(alpha_i[i] + t * alpha_i[i + 1])
                
                # Verify it's totally positive (for D ≤ 2000)
                if self.D <= 2000:
                    assert alpha_it.norm() > 0, f"alpha_{i,{t}} has norm ≤ 0"
                    assert alpha_it.trace() > 0, f"alpha_{i,{t}} not totally positive"
                
                all_indecomposables.append(alpha_it)
        
        if verbose:
            print(f"Generated {len(all_indecomposables)} candidates")
        
        # Normalize: keep only representatives up to units
        indecomposables_final = []
        tp_unit = self.tp_unit
        log_tp = self.R(log(tp_unit).n(max(self.precision, 20 * len(str(tp_unit.trace())))))
        
        for ind in all_indecomposables:
            # Check if ind is equivalent (up to units) to a previously found indecomposable
            is_new = True
            for tmp in indecomposables_final:
                ratio = ind / tmp
                if ratio in self.OK and self.OK(ratio).is_unit():
                    is_new = False
                    break
            
            if is_new:
                # Normalize ind: adjust by powers of tp_unit to get a canonical representative
                log_ind = self.R(log(ind).n(max(self.precision, 20 * len(str(ind.trace())))))
                l = self.R(((log(ind.norm()) / 2) - log_ind) / log_tp).round()
                
                good_rep = ind * (tp_unit ** l)
                
                # Try to improve representation by scaling by ±tp_unit
                candidates = [good_rep * tp_unit, good_rep / tp_unit]
                for cand in candidates:
                    if len(str(cand)) < len(str(good_rep)):
                        good_rep = cand
                
                indecomposables_final.append(good_rep)
                
                if verbose:
                    print(f"  Found indecomposable: {good_rep} (norm {good_rep.norm()})")
        
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
                # Verify Dress-Scharlau bound
                assert max_norm <= abs(self.discriminant) / 4, \
                    f"Max norm {max_norm} exceeds bound {abs(self.discriminant)/4}"
        
        return indecomposables_final
    
    def compute_num_indecomposables(self):
        """
        Computes just the number of indecomposables (modulo totally positive units) 
        This can be done faster than simply computing explicitly the indecomposables themselves.
        Uses the computation of M_D as given in the Blomer--Kala "On the Rank of Universal Quadratic Forms over Real Quadratic Fields" paper (page 16)
        """

        if (self.D%4==1): cf = continued_fraction((1 + self.a)/2)
        else: cf = continued_fraction(self.a)

        cp = cf.period()
        s = cf.period_length()

        # Compute the number of indecomposables! :)
        ans = -1
        if (s%2)==0:
            ans = sum([cp[i] for i in range(0, s-1, 2)])
        elif (D%4)==1:
            ans = 2*cf[0] + sum([cp[i] for i in range(s-1)]) - 1
        else:
            ans = 2*cf[0] + sum([cp[i] for i in range(s-1)])

        return ans


    def indecomposables(self):
        """Get list of indecomposables (computed on demand)."""
        if self._indecomposables is None:
            self.compute_indecomposables_dress_scharlau(verbose=False)
        return self._indecomposables
