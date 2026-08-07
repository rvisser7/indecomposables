"""
Explicit indecomposables in a biquadratic family due to Siu Hang Man.

Family:
    n >= 6,
    p = (2n-1)(2n+1),
    q = (2n-1)(2n+3),
    r = (2n+1)(2n+3),
    K = Q(sqrt(p), sqrt(q)).

This module mirrors the structure of simplest_cubic.py: generate explicit candidates,
then normalize up to multiplication by units.

        Note:
        The theorem statement as provided in project discussion contains one ambiguous line
        (a duplicated family involving eps_p and eps_r), and a range endpoint that appears
        inconsistent with the stated total count 10n-15. The implementation supports:
        - mode="as_stated": formulas exactly as provided,
        - mode="symmetric_guess": duplicated family replaced by natural q-analogue,
        - mode="count_formula_guess": symmetric guess and t-range in item (5) adjusted
            to match the stated total count 10n-15.
"""

from sage.all import Integer, PolynomialRing, QuadraticField


class BiquadraticField:
    """
    Explicit indecomposables for the Siu Hang Man family of biquadratic fields.

    The field is parameterized by n (not by arbitrary a,b), with
    K = Q(sqrt(p), sqrt(q)) and p,q,r as above.
    """

    def __init__(self, n, precision=100, validate_squarefree=True):
        self.n = Integer(n)
        self.precision = precision
        self.validate_squarefree = validate_squarefree

        if self.n < 6:
            raise ValueError("n must be at least 6")

        self.p = Integer((2 * self.n - 1) * (2 * self.n + 1))
        self.q = Integer((2 * self.n - 1) * (2 * self.n + 3))
        self.r = Integer((2 * self.n + 1) * (2 * self.n + 3))

        if self.validate_squarefree:
            for d in [self.p, self.q, self.r]:
                if not d.is_squarefree():
                    raise ValueError(f"Expected squarefree parameter, got d={d}")

        self.K = self._construct_field()
        self.OK = self.K.ring_of_integers()

        # Canonical square roots inside K.
        self.sqrt_p = self._sqrt_p
        self.sqrt_q = self._sqrt_q
        self.sqrt_r = self.sqrt_p * self.sqrt_q / Integer(2 * self.n - 1)

        self.discriminant = self.K.absolute_discriminant()
        self.regulator = self.K.regulator()
        self.class_number = self.K.class_number()

        self._indecomposables = None

    def _construct_field(self):
        """Construct K = Q(sqrt(p), sqrt(q)) as a quadratic tower."""
        Kp = QuadraticField(self.p, names="sp")
        self._sqrt_p = Kp.gen()

        Rp = PolynomialRing(Kp, "y")
        y = Rp.gen()
        K = Kp.extension(y**2 - self.q, names="sq")
        self._sqrt_q = K.gen()

        # Coerce sqrt(p) into the full field.
        self._sqrt_p = K(self._sqrt_p)
        return K

    @property
    def indecomposables(self):
        return self._indecomposables

    @staticmethod
    def _fundamental_unit_gt_one(d):
        """Return a fundamental unit eps of Q(sqrt(d)) with eps > 1."""
        F = QuadraticField(d, names="u")
        eps = F.unit_group().fundamental_units()[0]
        if eps.n(80) < 1:
            eps = eps**(-1)
        return F, eps

    def _embed_quadratic_element(self, element, sqrt_d_in_k):
        """Embed a + b*sqrt(d) from a quadratic subfield into K."""
        coeffs = element.list()
        if len(coeffs) != 2:
            raise ValueError(f"Expected quadratic element, got coefficients {coeffs}")
        return self.K(coeffs[0]) + self.K(coeffs[1]) * sqrt_d_in_k

    def _subfield_units_in_k(self):
        """Compute eps_p, eps_q, eps_r and embed them into K."""
        _, eps_p = self._fundamental_unit_gt_one(self.p)
        _, eps_q = self._fundamental_unit_gt_one(self.q)
        _, eps_r = self._fundamental_unit_gt_one(self.r)

        eps_p_k = self._embed_quadratic_element(eps_p, self.sqrt_p)
        eps_q_k = self._embed_quadratic_element(eps_q, self.sqrt_q)
        eps_r_k = self._embed_quadratic_element(eps_r, self.sqrt_r)

        return eps_p_k, eps_q_k, eps_r_k

    def _explicit_candidates(self, mode="count_formula_guess"):
        """
        Generate explicit candidates from the theorem statement.

        Args:
                        mode:
                                - "as_stated": implement the formulas exactly as provided.
                                - "symmetric_guess": replace duplicated family in item (5)
                                    by the natural q-analogue.
                                - "count_formula_guess": same as symmetric_guess, but use
                                    t = 4, ..., 2n-2 in item (5), matching total count 10n-15.
        """
        n = self.n
        K = self.K

        eps_p, eps_q, eps_r = self._subfield_units_in_k()

        mu = (
            K(Integer(2 * n + 3)) / 2
            + self.sqrt_p / 2
            + self.sqrt_q / 2
            + self.sqrt_r / 2
        )

        candidates = [
            K(1),
            (eps_p**(-1) + eps_r) / 2,
            mu,
        ]

        for t in range(2, 2 * n - 1):
            candidates.append(K(1) + eps_p + t * (mu - K(1)))
            candidates.append((K(1) + eps_p * eps_r) / 2 + eps_p + t * (mu - K(1)))

        t_upper = 2 * n
        if mode == "count_formula_guess":
            t_upper = 2 * n - 1

        for t in range(4, t_upper):
            candidates.append(K(1) + eps_q**(-1) + t * (mu - K(1)))
            if mode == "as_stated":
                candidates.append((K(1) + eps_p * eps_r) / 2 + eps_p + t * (mu - K(1)))
            elif mode in ["symmetric_guess", "count_formula_guess"]:
                candidates.append((K(1) + eps_q * eps_r) / 2 + eps_q + t * (mu - K(1)))
            else:
                raise ValueError(f"Unsupported mode={mode}")

        for t in range(2, 2 * n):
            candidates.append((eps_r**(-1) + eps_p) / 2 + t * (mu - K(2)))

        return candidates

    def compute_indecomposables_man(
        self,
        mode="count_formula_guess",
        verbose=True,
        check_theorem_count=False,
        require_integral=False,
        require_totally_positive=False,
        reduce_mod_units=False,
    ):
        """
        Compute theorem representatives for indecomposables.

        By default this returns the explicit representatives listed by the theorem.
        Optional flags can enforce integrality / total positivity and reduction modulo units.
        """
        if self._indecomposables is not None:
            return self._indecomposables

        candidates = self._explicit_candidates(mode=mode)

        if verbose:
            print(f"Computing Man-family indecomposables for n={self.n}, mode={mode}")
            print(f"Generated {len(candidates)} explicit representatives")

        filtered = []
        for c in candidates:
            cc = self.K(c)
            if require_integral:
                try:
                    cc = self.K(self.OK(cc))
                except Exception:
                    continue
            if require_totally_positive and not cc.is_totally_positive():
                continue
            filtered.append(cc)

        indecomposables_final = []
        if reduce_mod_units:
            for cand in filtered:
                is_new = True
                for prev in indecomposables_final:
                    ratio = cand / prev
                    if ratio in self.OK and self.OK(ratio).is_unit():
                        is_new = False
                        break
                if is_new:
                    indecomposables_final.append(cand)
        else:
            # Keep theorem representatives in their listed order, removing exact duplicates only.
            seen = set()
            for cand in filtered:
                key = str(cand)
                if key not in seen:
                    seen.add(key)
                    indecomposables_final.append(cand)

        indecomposables_final.sort(
            key=lambda x: (self.K(x).norm(), self.K(x).trace(), len(str(x)))
        )

        self._indecomposables = indecomposables_final

        if check_theorem_count:
            expected = self.theorem_expected_count()
            if len(indecomposables_final) != expected:
                raise ValueError(
                    f"Count mismatch: got {len(indecomposables_final)}, expected {expected}. "
                    "Try mode='symmetric_guess' if using a corrected theorem variant."
                )

        if verbose:
            print(f"Total representatives returned: {len(indecomposables_final)}")

        return indecomposables_final

    def theorem_expected_count(self):
        """Expected count from the theorem statement: 10n - 15."""
        return Integer(10 * self.n - 15)

    def to_data_row(self):
        """Serialize result in pipe-delimited row (similar to other scripts)."""
        if self._indecomposables is None:
            raise ValueError("Indecomposables not computed yet")

        ks = self._indecomposables
        min_norm = min(self.K(k).norm() for k in ks) if ks else 0
        max_norm = max(self.K(k).norm() for k in ks) if ks else 0

        return (
            f"biquad_man_n={self.n}|{self.discriminant}|{self.regulator}|{self.class_number}|"
            f"{len(ks)}|{min_norm}|{max_norm}|[{', '.join(str(x) for x in ks)}]"
        )


if __name__ == "__main__":
    bq = BiquadraticField(6)
    vals = bq.compute_indecomposables_man(mode="count_formula_guess", verbose=True)
    print(f"Expected theorem count: {bq.theorem_expected_count()}")
    print(f"Computed count: {len(vals)}")
