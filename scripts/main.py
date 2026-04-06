# Main Python file for computing indecomposables in totally real number fields

# Import Sage modules
from sage.all import (
    Set, ZZ, RR, pi, euler_phi, CyclotomicField, gap, RealField, sqrt, prod,
    QQ, NumberField, PolynomialRing, latex, pari, cached_function, Permutation,
    sign, Matrix, GF, log, proof, QuadraticField, AA)

from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens
import itertools

# Import RealQuadraticField
from real_quadratic import RealQuadraticField as RQ
# Import SimplestCubicField
from simplest_cubic import SimplestCubicField as SCF

from logging_utils import setup_logger

logger = setup_logger(
    log_file=None,
    verbose=False
)

ZZx = PolynomialRing(ZZ, names='x')

# Sage's .is_totally_positive() uses AA (algebraic real field), so should be safe to use
# e.g. see: https://github.com/sagemath/sage/blob/31a24ce25741e1610aa90924ce637c018ee6d87d/src/sage/rings/number_field/number_field_element.pyx#L2139

# At the moment, we're only really supporting totally real fields
# Should we maybe call this class "TotallyRealField" or "TotallyRealFieldData" instead??
class NumberFieldData:
    """
    Class for computing and storing indecomposable elements and sails in totally real number fields.
    
    An element alpha in a number field K is totally positive if σ(α) > 0 for all real embeddings σ.
    An element alpha in O_K is indecomposable if it cannot be written as a sum of two other 
    totally positive integral elements in O_K. 
    """
    
    def __init__(self, coeffs, lmfdb_label=None, lmfdb_index=None, metadata=None):
        """
        Initialize NumberFieldData.
        
        Args:
            lmfdb_label: LMFDB label for the field (e.g., "3.3.49.1")
            lmfdb_index: Just the last part of the LMFDB label (i.e. "1")
            coeffs:  Coefficients of a defining polynomial for the field
            metadata: Dict with field metadata (degree, discriminant, regulator, class_number, etc.)
        """
        logger.debug("Initialisaing new number field class")
                    
        self.lmfdb_label = lmfdb_label
        self.lmfdb_index = lmfdb_index
        self.coeffs = coeffs
        self.metadata = metadata or {}

        # Construct the field
        if self.coeffs is not None:
            self.K = NumberField(ZZx(coeffs), names='a')
        
        # Computed data (populated by compute_* methods)
        self._degree = None
        self._discriminant = None
        self._regulator = None
        self._class_number = None
        self._fundamental_units = None
        self._totally_positive_units = None
        self._unit_representatives = None
        self._indecomposables = None
        self._indecomposables_logs = None
        self._unit_signature_rank = None

        # A dictionary which converts a signature into a unit contained in that signature (if it exists)
        self._unit_in_signature = {}
        self._tp_unit_index = None     # kept around for historical reasons
        
        # Configuration (should ideally phase these out...)
        self.exit_for_nonunit = False  # Exit early if non-unit indecomposable found
        self.check_asserts = False     # Run expensive assertions
        self.max_norm_bound = None     # Override default norm bound
        self.precision = 100           # RealField precision
    
    @property
    def discriminant(self):
        """Get discriminant (computed if not cached)."""
        if self._discriminant is None and self.K is not None:
            self._discriminant = self.K.discriminant()
        return self._discriminant or self.metadata.get('discriminant')
    
    @property
    def regulator(self):
        """Get regulator (computed if not cached)."""
        if self._regulator is None and self.K is not None:
            self._regulator = self.K.regulator()
        return self._regulator or self.metadata.get('regulator')
    
    @property
    def class_number(self):
        """Get class number (computed if not cached)."""
        if self._class_number is None and self.K is not None:
            self._class_number = self.K.class_number()
        return self._class_number or self.metadata.get('class_number')
    
    @property
    def degree(self):
        """Get degree of the number field."""
        if self.K is not None:
            return self.K.degree()
        return self.metadata.get('degree')
    
    @property
    def indecomposables(self):
        """Get list of indecomposables (up to multiplication by totally positive units)."""
        return self._indecomposables
    
    @property
    def fundamental_units(self):
        """Get fundamental units of the field."""
        if self._fundamental_units is None and self.K is not None:
            self._compute_fundamental_units()
        return self._fundamental_units
    
    @property
    def totally_positive_units(self):
        """Get basis for totally positive units."""
        if self._totally_positive_units is None and self.K is not None:
            self._compute_totally_positive_units()
        return self._totally_positive_units
    
    @property
    def unit_representatives(self):
        """Get representatives of all unit classes modulo squares."""
        if self._unit_representatives is None and self.K is not None:
            self._compute_unit_representatives()
        return self._unit_representatives


    def _signature(self, x):
        """ Computes the signature of an element x in this field K """

        embeds = (self.K(x)).embeddings(AA)
        return tuple([e.sign() for e in embeds])
        
    
    def _compute_fundamental_units(self):
        """Compute a set of fundamental units for the field."""
        logger.info("Computing set of fundamental units")

        if self.K is None:
            raise ValueError("Number field K must be set first")
        
        UK = self.K.unit_group()
        fun_units = UK.fundamental_units()
        
        # Verify they're units
        #assert all(u.norm() in [1, -1] for u in fun_units)
        
        # Normalize so all positive w.r.t. first embedding
        normalized_units = []
        for u in fun_units:
            u = self.K(u)
            embeddings = u.complex_embeddings(self.precision)
            if embeddings[0] < 0:
                u = -u
            normalized_units.append(u)
            # Verify positivity under first embedding
            assert self.K(u).complex_embeddings(self.precision)[0] > 0
        
        self._fundamental_units = normalized_units
        
    
    def _compute_totally_positive_units(self):
        """
        Compute a basis for the group of totally positive units.
        Also computes the unit signature rank.
        
        Uses the approach from https://mathoverflow.net/questions/483102/
        to solve for the kernel of a GF(2) matrix.
        """
        logger.debug("Computing totally positive units...")

        if self._fundamental_units is None:
            self._compute_fundamental_units()
        
        d = self.degree
        fun_units = self._fundamental_units
        
        # Set up matrix over GF(2): M[i][j] = 1 if σ_i(u_j) < 0
        F2 = GF(2)
        M = [[F2(0) for _ in range(d-1)] for _ in range(d-1)]
        
        for i in range(d-1):
            for j in range(d-1):
                embeddings = self.K(fun_units[j]).complex_embeddings(self.precision)
                if embeddings[i+1] < 0:  # Skip first embedding
                    M[i][j] = F2(1)
        
        M = Matrix(GF(2), M)
        kernel = M.right_kernel()
        kernel_basis = kernel.basis()
        
        # Build totally positive units from kernel
        utp = []
        utp_exponents = []
        
        # First: from kernel basis
        for kk in kernel_basis:
            new_unit = prod([fun_units[i]**(ZZ(kk[i])) for i in range(d-1)])
            utp.append(self.K(new_unit))
            utp_exponents.append([ZZ(t) for t in kk])
        
        # Then: add u_i^2 if needed to reach rank d-1
        for i in range(d-1):
            v_to_add = [ZZ(0)] * (d-1)
            v_to_add[i] = ZZ(2)
            
            to_add = True
            if len(utp) > 0:
                Mtmp = Matrix(ZZ, len(utp), d-1, utp_exponents)
                try:
                    Mtmp.solve_left(Matrix(ZZ, 1, d-1, v_to_add))
                    to_add = False
                except ValueError:
                    to_add = True
            
            if to_add:
                new_unit = fun_units[i]**2
                utp.append(self.K(new_unit))
                utp_exponents.append(v_to_add)
        
        assert len(utp) == d-1, f"Expected {d-1} totally positive units, got {len(utp)}"
        
        # Verify all are totally positive units
        for u in utp:
            assert u.norm() == 1, f"Unit {u} has norm {u.norm()}"
            assert u.is_totally_positive()
        
        self._totally_positive_units = utp
        self._tp_unit_index = 2**(d-1 - len(kernel_basis))
        self._unit_signature_rank = len(kernel_basis)
    
    def _compute_unit_representatives(self):
        """
        Compute representatives of the unit group U_K modulo squares of untis (U_K^2).
        Returns a set of size 2^(r-1)
        """
        logger.info("Computing unit representatives...")

        d = self.degree
        fun_units = self.fundamental_units
        
        representatives = []
        for i in range(pow(2, d)):
            tmp = self.K.ring_of_integers()(1)
            
            # Handle ±1
            if (i % 2) == 0:
                tmp *= self.K.ring_of_integers()(-1)
            
            # Handle fundamental units
            for j in range(1, d):
                if ((i // (2**j)) % 2) == 1:
                    tmp *= fun_units[j-1]
            
            # Verify it's a unit
            assert tmp.norm() in [1, -1]
            assert tmp in self.K.ring_of_integers()
            
            representatives.append(tmp)
        
        self._unit_representatives = representatives
    
    def _compute_max_norm_bound(self):
        """
        Compute bound for maximal norm of indecomposable based on field degree.
        
        Uses Kala-Yatsyna bound and field-specific optimizations.
        Reference: Theorem 5 of Kala, Vítězslav; Yatsyna, Pavlo; 
        On Kitaoka's conjecture and lifting problem for universal quadratic forms.
        """
        if self.max_norm_bound is not None:
            return self.max_norm_bound
        
        disc = self.discriminant
        d = self.degree
        
        if d == 1:
            return 1
        elif d == 2:
            return abs(disc) // 4 
        elif d == 3:
            if disc <= 1000:
                return abs(disc)
            elif disc <= 10000:
                return abs(disc) // 5
            elif disc <= 50000:
                return abs(disc) // 10
            else:
                return abs(disc) // 20
        elif d == 4:
            return min(abs(disc), 100000)
        elif d == 5:
            return abs(disc) // 10
        elif d == 6:
            return 10000
        else:
            return 10000
    
    def _ideal_gens_of_norm(self, n, algorithm="new"):
        """
        Returns generators for principal ideals of norm n

        Here, if algorithm is "old", then uses the built-in K.ideals_of_bdd_norm
        Otherwise, if algorithm is "new", then uses the Pavlo-suggested bdd_norm_pr_ideal_gens function
        """

        if algorithm=="old":
            ans = []
            idls = self._ideals_of_bdd_norm[n]
            for I in idls:
                gens = I.gens_reduced()
                if len(gens) != 1:
                    continue
                ans.append(gens[0])
            return ans

        if algorithm=="new":
            return bdd_norm_pr_ideal_gens(self.K, [n])
    


    def compute_indecomposables_brute_force(self, verbose=True, algorithm="new", GRH=False):
        """
        Main brute force algorithm: compute all indecomposables up to totally positive units.
        
        Uses brute force enumeration over ideals of bounded norm. For each principal
        ideal with a totally positive generator, checks if it can be decomposed as
        a sum of previously found indecomposables.

        Uses the ideal bound brute force method (based on Kala-Yatsyna bound) to enumerate
        all indecomposables up to multiplication by totally positive units.
        
        Args:
            verbose: Print progress information
        
        Returns:
            List of indecomposables, sorted by norm then trace
        """
        if verbose:
            print("Starting brute force indecomposables computation for field K")

        if self.K is None:
            raise ValueError("Number field K must be set first")
        
        # Assume GRH (to make computations faster)
        if GRH:
            proof.number_field(False)
        
        # Prepare units
        if self._fundamental_units is None:
            self._compute_fundamental_units()
        if self._totally_positive_units is None:
            self._compute_totally_positive_units()
        if self._unit_representatives is None:
            self._compute_unit_representatives()
        
        d = self.degree
        OK = self.K.ring_of_integers()
        disc = self.discriminant
        utp = self.totally_positive_units
        unit_reps = self.unit_representatives
        R = RealField(self.precision)
        
        max_norm_bound = self._compute_max_norm_bound()
        if verbose:
            print(f"Discriminant: {disc}")
            print(f"Max norm bound for search: {max_norm_bound}")
        
        # Build log matrix for unit group
        M = [[R(0) for _ in range(d-1)] for _ in range(d)]
        for i in range(d-1):
            prec = max(self.precision, 20 * len(str(max(abs(c) for c in utp[i].minpoly().coefficients()))))
            for j in range(d):
                emb = self.K(utp[i]).complex_embeddings(prec)[j]
                M[j][i] = R(emb.log())
        
        # Storage for indecomposables and their logs
        all_indecomposables = []
        all_indecomposables_logs = []
        
        # Compute indecomposables by checking ideals of bounded norm
        max_k_witness = -1
        num_checked = 0
        
        for nrm in range(1, max_norm_bound + 1):
            if verbose:
                if ((nrm%10000)==0) or ((nrm < 100000) and ((nrm%1000)==0)) or ((nrm < 1000) and ((nrm%100)==0)):
                    percent = 100.0 * nrm / max_norm_bound
                    print(f"Checking ideals of norm {nrm} ({percent:.1f}% to max norm bound), largest k-witness so far: {max_k_witness})")
            
            # Get ideal generators of this norm
            ideal_gens_of_norm = self._ideal_gens_of_norm(nrm, algorithm=algorithm)
                        
            for gen in ideal_gens_of_norm:
                num_checked += 1
                
                # Find totally positive generator
                tp_gen = self._find_totally_positive_generator(gen, unit_reps)
                if tp_gen is None:
                    continue
                
                # Adjust to minimize coefficients
                tp_gen = self._normalize_totally_positive_gen(tp_gen, utp, M, nrm)
                
                # Compute logs of tp_gen
                prec = max(self.precision, 20 * len(str(max(abs(c) for c in tp_gen.minpoly().coefficients()))))
                log_tp_gen = [R(tp_gen.complex_embeddings(prec)[i].log()) for i in range(d)]
                
                # Check if decomposable
                is_indecomp = self._is_indecomposable_unit_algorithm(
                    tp_gen, log_tp_gen, all_indecomposables, all_indecomposables_logs,
                    utp, M, nrm, R
                )
                
                if is_indecomp:
                    all_indecomposables.append(tp_gen)
                    all_indecomposables_logs.append(log_tp_gen)
                    if verbose:
                        print(f"  Found indecomposable: {tp_gen} (norm {self.K(tp_gen).norm()})")
                    
                    if self.exit_for_nonunit and nrm > 1:
                        break
            
            if self.exit_for_nonunit and len(all_indecomposables) > 1:
                break
        
        # Sort by norm, then trace, then string representation
        all_indecomposables.sort(key=lambda x: (self.K(x).norm(), self.K(x).trace(), len(str(x))))
        
        self._indecomposables = all_indecomposables
        self._indecomposables_logs = all_indecomposables_logs
        
        if verbose:
            print(f"\nTotal indecomposables found: {len(all_indecomposables)}")
        
        return all_indecomposables
    
    def compute_indecomposables(self, verbose=True):
        """
        Compute indecomposables using whichever is the fastest algorithm available.
        
        For real quadratic fields (degree 2), uses the Dress-Scharlau continued fraction
        method (much faster than the general brute force approach).
        
        For simplest cubic fields (degree 3) where Z[rho] has index 1 or 3 in the maximal order,
        uses the Kala-Tinkovaa and Gil-Munoz--Tinkova classification (much faster than brute force).

        For certain biquadratics, can use the classification by Siu Hang Man.
        
        For other fields, falls back to the general brute force method.
        
        Args:
            verbose: Print progress information
        
        Returns:
            List of indecomposables (up to multiplication by totally positive units)
        
        """
            
        if self.K is None:
            raise ValueError("Number field K must be set first")
        
        # Use Dress-Scharlau for real quadratic fields
        if self.degree == 2:
            if verbose:
                print("Using Dress-Scharlau method for real quadratic field")
            
            # Extract D from the quadratic field
            # For a quadratic field defined by X^2 - D or similar
            # Try to create RealQuadraticField

            D = self.discriminant.squarefree_part()
            rq = RQ(D, precision=self.precision)
            self._indecomposables = rq.compute_indecomposables_dress_scharlau(verbose=verbose)
            return self._indecomposables
                    
        # Use Kala-Tinkova for simplest cubic fields
        if self.degree == 3:
            if verbose:
                print("Checking if field is a simplest cubic field...")
            
            try:
                # Check if this is a simplest cubic field of the form x^3 - n*x^2 - (n+3)*x - 1
                poly = self.K.defining_polynomial()
                if poly.degree() == 3:
                    # Check if it matches the simplest cubic form
                    # x^3 - n*x^2 - (n+3)*x - 1
                    a3, a2, a1, a0 = poly.list()
                    if a3 == 1 and a0 == -1 and a1 == -(a2 + 3):
                        n = -a2  # Since -n*x^2, so a2 = -n
                        if n in ZZ and n >= -1:
                            scf = SCF(n, precision=self.precision)
                            if scf.can_use_classification:
                                if verbose:
                                    print(f"Using Kala-Tinkova classification for simplest cubic field with n = {n}")
                                self._indecomposables = scf.compute_indecomposables_kala_tinkova(verbose=verbose)
                                return self._indecomposables
                            else:
                                if verbose:
                                    print(f"Simplest cubic field with n = {n} has index {scf.index}, cannot use classification")
            except Exception as e:
                if verbose:
                    print(f"Could not use Kala-Tinková: {e}")
                    print("Falling back to brute force method")
        
        # Fall back to general brute force
        if verbose:
            print("Using general brute force method")
        return self.compute_indecomposables_brute_force(verbose=verbose)
    
    def _find_totally_positive_generator(self, gen, unit_reps):
        """Find a totally positive generator by multiplying by units."""
        for u in unit_reps:
            tmp = gen * u
            embeddings = self.K(tmp).complex_embeddings(self.precision)
            if all(e > 0 for e in embeddings):
                return tmp
        return None
    
    def _normalize_totally_positive_gen(self, tp_gen, utp, M, nrm):
        """
        Normalize tp_gen by adjusting with totally positive units to minimize coefficients.
        
        Solves: log(nrm)/d - log(tp_gen) = M * l
        to find optimal l, then returns tp_gen * prod(utp[i]^l[i])
        """
        d = self.degree
        R = RealField(self.precision)
        prec = max(self.precision, 20 * len(str(max(abs(c) for c in tp_gen.minpoly().coefficients()))))
        
        M_red = Matrix(R, M[:-1])  # Remove last row
        log_tp_gen = [R(tp_gen.complex_embeddings(prec)[i].log()) for i in range(d)]
        V_right = Matrix(R, d-1, 1, [R(log(nrm)/d) - log_tp_gen[i] for i in range(d-1)])
        
        try:
            sol = M_red.solve_right(V_right)
            l = [tmp[0].round() for tmp in sol]
            tp_gen = tp_gen * prod([utp[i]**l[i] for i in range(d-1)])
        except (ZeroDivisionError, ValueError):
            pass  # Return original if solve fails
        
        return tp_gen
    
    def _is_indecomposable_unit_algorithm(self, tp_gen, log_tp_gen, all_indecomposables, 
                               all_indecomposables_logs, utp, M, nrm, R):
        """
        Check if tp_gen is indecomposable by testing against all previously found indecomposables.
        and for each one, solve an integer feasbility problem using the unit log-lattice embedding
        
        For each indecomposable ind, tries to find k ∈ Z^(d-1) such that
        tp_gen = ind * prod(utp[i]^k[i]) + something_totally_positive
        """
        d = self.degree
        prec = max(self.precision, 20 * len(str(max(abs(c) for c in tp_gen.minpoly().coefficients()))))
        
        for ind, log_ind in zip(all_indecomposables, all_indecomposables_logs):
            # Use Cramer's rule on various pairs of embeddings to bound k
            k_lower = [1000000 for _ in range(d-1)]
            k_upper = [-1000000 for _ in range(d-1)]
            
            for case in range(d):
                # Remove case-th row from M
                M_tmp = Matrix(R, M[:case] + M[case+1:])
                V_right = Matrix(R, d-1, 1, 
                                [log_tp_gen[j] - log_ind[j] 
                                 for j in list(range(case)) + list(range(case+1, d))])
                
                sol = M_tmp.solve_right(V_right)
                for j in range(d-1):
                    k_lower[j] = min(k_lower[j], sol[j][0].floor())
                    k_upper[j] = max(k_upper[j], sol[j][0].ceil())
            
            # Check all possible k in bounded range
            iter_ranges = [range(k_lower[j]-0, k_upper[j]+1) for j in range(d-1)]
            for k in itertools.product(*iter_ranges):
                ind_shifted = ind * prod([utp[i]**k[i] for i in range(d-1)])
                diff = tp_gen - ind_shifted
                
                embeddings = self.K(diff).complex_embeddings(prec)
                if all(e > 0 for e in embeddings):
                    # Found decomposition
                    return False
        
        # No decomposition found
        return nrm < 2**d or True  # Force True for now; can adjust logic
    

    def compute_sails(self, signature, verbose=True, GRH=False):
        """
        Algorithm to compute the sails of K
        """

        # TODO
        pass


    def to_data_row(self):
        """
        Format as a pipe-delimited data row for output file.
        
        Returns string with columns: lmfdb_label|discriminant|regulator|class_number|etc.
        """
        if self._indecomposables is None:
            raise ValueError("Must compute indecomposables first")
        
        d = self.degree
        indecomp_list = self._indecomposables
        
        # Get min norm of non-unit indecomposables
        min_norm = 1000000
        if len(indecomp_list) > 1:
            min_norm = self.K(indecomp_list[1]).norm()
        
        max_norm = self.K(indecomp_list[-1]).norm() if indecomp_list else 1
        
        # Build row
        row = f"{self.lmfdb_label}|"
        row += f"{self.discriminant}|"
        row += f"{self.regulator}|"
        row += f"{self.class_number}|"
        row += f"{len(indecomp_list)}|"
        row += f"{min_norm}|"
        row += f"{max_norm}|"
        row += "["
        row += ", ".join(str(ind) for ind in indecomp_list)
        row += "]"
        
        return row


    def from_data_row(self):
        """
        Populate properties of this field from data given in a pipe-delimited data row
        (probably shouldn't be used too often)
        """

        # TODO
        pass


# ============================================================================
# Utility Functions
# ============================================================================

def create_field_from_polynomial(poly, var='a'):
    """
    Create a NumberFieldData from a polynomial.
    
    Args:
        poly: A polynomial in Sage/Python (or a string like 'x^3 - 2')
        var: Variable name to use (default 'a' for the primitive element)
    
    Returns:
        NumberFieldData object with K set to the number field
    """
    from sage.all import PolynomialRing, QQ as QQ_sage
    
    if isinstance(poly, str):
        R = PolynomialRing(QQ_sage, 'x')
        x = R.gen()
        # Parse string polynomial (simple eval, consider using sage.calculus.calculus for safety)
        poly = eval(poly.replace('x', 'x'))
    
    K = poly.parent().base_ring().extension(poly, names=var)
    return NumberFieldData(field=K, metadata={'degree': K.degree()})


def process_degree_file(filename):
    """
    Process a data file to create NumberFieldData objects for each field.
    
    Assumes format: header line, type line, blank line, then data rows.
    
    Args:
        filename: Path to data file (e.g., "data/deg3.txt")
    
    Returns:
        List of NumberFieldData objects
    """
    fields = []
    
    with open(filename, 'r') as f:
        # Skip header and type line
        header = f.readline().strip().split('|')
        types = f.readline().strip().split('|')
        f.readline()  # Skip blank line
        
        # Process data rows
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            data = line.split('|')
            record = {h: d for h, d in zip(header, data)}
            
            # Create NumberFieldData
            nfd = NumberFieldData(
                label=record.get('lmfdb_label'),
                metadata={
                    'discriminant': int(record.get('discriminant', 0)),
                    'regulator': float(record.get('regulator', 0.0)),
                    'class_number': int(record.get('class_number', 0)),
                    'degree': 0  # Will be set when field is constructed
                }
            )
            fields.append(nfd)
    
    return fields


def batch_compute_indecomposables(fields, degree=3, verbose=True, output_file=None):
    """
    Compute indecomposables for a batch of fields.
    
    Args:
        fields: List of NumberFieldData objects
        degree: Degree of fields (used for configuration)
        verbose: Print progress
        output_file: Optional file to write output rows
    
    Returns:
        List of NumberFieldData objects with computed indecomposables
    """
    results = []
    
    for i, field_data in enumerate(fields):
        if verbose:
            print(f"\n[{i+1}/{len(fields)}] Processing {field_data.label}...")
        
        try:
            indecomp = field_data.compute_indecomposables(verbose=verbose)
            results.append(field_data)
            
            if output_file:
                row = field_data.to_data_row()
                with open(output_file, 'a') as f:
                    f.write(row + '\n')
        
        except Exception as e:
            if verbose:
                print(f"  Error: {e}")
    
    return results


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    # Example 1: Create a field from a polynomial and compute indecomposables
    print("=" * 60)
    print("Example 1: Cubic field x^3 - 2")
    print("=" * 60)
    
    from sage.all import PolynomialRing, QQ
    
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    poly = x**3 - 2
    
    # Create the number field
    K = QQ.extension(poly, names='a')
    print(f"Number field: {K}")
    print(f"Discriminant: {K.discriminant()}")
    print(f"Degree: {K.degree()}")
    
    # Create NumberFieldData object
    nfd = NumberFieldData(label="3.3.49.1", field=K)
    
    # Compute indecomposables
    try:
        indecomp = nfd.compute_indecomposables(verbose=True)
        print(f"\nFound {len(indecomp)} indecomposable(s):")
        for ind in indecomp:
            print(f"  {ind} (norm: {K(ind).norm()})")
    except Exception as e:
        print(f"Error during computation: {e}")
        import traceback
        traceback.print_exc()


