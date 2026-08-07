"""
Canonical representatives for totally positive elements modulo totally positive units.

The convention implemented here is:

    Among the U^+-orbit of x, take the elements of minimal trace; among those,
    take the one whose coordinate vector with respect to a fixed integral basis
    is lexicographically least.

This is well defined: for y totally positive, 0 < sigma_i(y) < Tr(y) for every i,
so each orbit contains only finitely many elements of bounded trace and the
minimum is attained.  The action of U^+ is free (a totally real field has only
+-1 as roots of unity, and -1 is not totally positive), so the minimiser is a
single element once ties are broken.

Design note on exactness
------------------------
Floating point is used *only to generate a candidate set*.  The winner is chosen
by comparing traces exactly in the order.  So precision can only ever affect
completeness of the candidate set, never the correctness of the comparison.
``canonical`` escalates precision until the answer is stable, so a
precision-starved run costs time rather than giving a wrong answer.

Signs (which determine total positivity and the unit signature map) are computed
in AA, i.e. exactly -- never from floating point embeddings.
"""

from sage.all import (
    AA,
    GF,
    Matrix,
    QQ,
    RealField,
    ZZ,
    identity_matrix,
    vector,
)

__all__ = [
    "UnitOrbitNormalizer",
    "totally_positive_unit_basis",
]


# ---------------------------------------------------------------------------
# Totally positive units
# ---------------------------------------------------------------------------

def totally_positive_unit_basis(K):
    """
    Return a Z-basis for the group U^+ of totally positive units of K.

    K is assumed totally real, so U has rank n-1 and torsion {+-1}, and U^+ is
    free of rank n-1.

    An element eps = prod u_i^{a_i} lies in +-U^+ exactly when its signature
    vector is constant.  So we take the kernel of

        Z^r --> F_2^n / <(1,...,1)>,   a |--> signature(prod u_i^{a_i})

    and fix the sign of each resulting generator.

    INPUT:

    - ``K`` -- a totally real number field

    OUTPUT: list of ``n-1`` totally positive units generating U^+

    EXAMPLES::

        sage: from normalize import totally_positive_unit_basis
        sage: K.<a> = NumberField(x^2 - 5)
        sage: B = totally_positive_unit_basis(K)
        sage: len(B) == 1 and all(u.is_totally_positive() for u in B)
        True
    """
    n = K.degree()
    embs = K.embeddings(AA)
    if len(embs) != n:
        raise ValueError(f"{K} is not totally real")

    units = list(K.unit_group().fundamental_units())
    r = len(units)
    if r != n - 1:
        raise ValueError(f"expected unit rank {n - 1}, got {r}")

    def sig(u):
        # exact signs via AA
        return [ZZ(0) if e(u) > 0 else ZZ(1) for e in embs]

    # Quotient F_2^n by the all-ones vector: v |--> (v_0+v_1, ..., v_0+v_{n-1}).
    rows = []
    for u in units:
        s = sig(u)
        rows.append([s[0] + s[i] for i in range(1, n)])
    M = Matrix(GF(2), r, n - 1, rows)

    # Lattice L = {a in Z^r : a*M = 0 mod 2}, spanned by lifted kernel + 2*Z^r.
    ker = [vector(ZZ, [ZZ(c) for c in v]) for v in M.left_kernel().basis()]
    gens = [list(v) for v in ker] + [list(v) for v in 2 * identity_matrix(ZZ, r)]
    L = Matrix(ZZ, gens).hermite_form(include_zero_rows=False)

    basis = []
    for row in L.rows():
        eps = K.one()
        for u, e in zip(units, row):
            eps *= u ** ZZ(e)
        if not eps.is_totally_positive():
            eps = -eps
        if not eps.is_totally_positive():
            raise RuntimeError("kernel computation produced a non-constant signature")
        basis.append(eps)

    if len(basis) != n - 1:
        raise RuntimeError(f"expected rank {n - 1} for U^+, got {len(basis)}")
    return basis


# ---------------------------------------------------------------------------
# Lattice enumeration (Fincke-Pohst, shifted)
# ---------------------------------------------------------------------------

def _cholesky_q(G, r):
    """Cohen Alg. 2.7.6: Q(z) = sum_i q[i][i] * (z_i + sum_{j>i} q[i][j] z_j)^2."""
    q = [[G[i][j] for j in range(r)] for i in range(r)]
    for i in range(r):
        if q[i][i] <= 0:
            raise ValueError("Gram matrix not positive definite at working precision")
        for j in range(i + 1, r):
            q[j][i] = q[i][j]
            q[i][j] = q[i][j] / q[i][i]
        for k in range(i + 1, r):
            for ll in range(k, r):
                q[k][ll] = q[k][ll] - q[k][i] * q[i][ll]
    return q


def _enumerate_close(G, y, rho2, r, cap):
    """
    All a in Z^r with (a - y) G (a - y)^T <= rho2.

    Returns a *superset* is acceptable to callers: integer ranges are widened
    outward so that floating point can only ever add candidates, not drop them.
    """
    q = _cholesky_q(G, r)
    out = []
    a = [ZZ(0)] * r
    slack = G.base_ring()(2) ** (-G.base_ring().precision() // 2)

    def rec(i, T):
        if len(out) > cap:
            raise OverflowError(
                f"more than {cap} candidate units; the log lattice is probably "
                "badly reduced or the working precision is too low")
        if T < -slack:
            return
        U = sum(q[i][j] * (a[j] - y[j]) for j in range(i + 1, r))
        rad2 = T / q[i][i]
        if rad2 < 0:
            rad2 = rad2.parent()(0)
        rad = rad2.sqrt() + slack
        centre = y[i] - U
        lo = (centre - rad).floor()
        hi = (centre + rad).ceil()
        for ai in range(lo, hi + 1):
            a[i] = ZZ(ai)
            term = q[i][i] * (a[i] - y[i] + U) ** 2
            if i == 0:
                if term <= T + slack:
                    out.append(tuple(a))
            else:
                rec(i - 1, T - term)
        a[i] = ZZ(0)

    rec(r - 1, rho2)
    return out


# ---------------------------------------------------------------------------
# Normalizer
# ---------------------------------------------------------------------------

class UnitOrbitNormalizer:
    """
    Canonical representatives for the action of U^+ on totally positive elements.

    One instance per field.  Constructing it computes the unit group and the
    LLL-reduced log lattice, which is the expensive part, so instances should be
    cached per field and reused across all indecomposables of that field.

    INPUT:

    - ``K`` -- a totally real number field
    - ``basis`` -- optional integral basis (default: basis of the maximal order).
      This fixes the lexicographic tie-break and the stored coordinate vectors,
      so it is part of the published schema: use LMFDB's ``zk``.
    - ``unit_basis`` -- optional precomputed Z-basis of U^+
    - ``prec`` -- initial working precision in bits (default 128)

    EXAMPLES::

        sage: from normalize import UnitOrbitNormalizer
        sage: K.<a> = NumberField(x^2 - 5)
        sage: N = UnitOrbitNormalizer(K)
        sage: eps = N.unit_basis[0]
        sage: x = K(11)
        sage: N.canonical(eps^3 * x) == N.canonical(x)
        True
        sage: N.canonical(eps^7) == K.one()
        True
    """

    MAX_CANDIDATES = 200000

    def __init__(self, K, basis=None, unit_basis=None, prec=128):
        self.K = K
        self.n = K.degree()
        self.prec = prec

        self.embeddings_exact = K.embeddings(AA)
        if len(self.embeddings_exact) != self.n:
            raise ValueError(f"{K} is not totally real")

        self.basis = list(basis) if basis is not None else list(K.maximal_order().basis())
        if len(self.basis) != self.n:
            raise ValueError("integral basis has the wrong length")
        self._basis_inv = Matrix(QQ, [K(b).vector() for b in self.basis]).inverse()

        self.unit_basis = (list(unit_basis) if unit_basis is not None
                           else totally_positive_unit_basis(K))
        self.rank = len(self.unit_basis)

    # -- coordinates -------------------------------------------------------

    def coordinates(self, x):
        """Coordinate vector of ``x`` with respect to the chosen integral basis."""
        c = self.K(x).vector() * self._basis_inv
        if not all(t in ZZ for t in c):
            raise ValueError(f"{x} is not integral with respect to the chosen basis")
        return vector(ZZ, [ZZ(t) for t in c])

    def sort_key(self, x):
        """Total order on elements: norm, then trace, then coordinates."""
        x = self.K(x)
        return (x.norm(), x.trace(), tuple(self.coordinates(x)))

    # -- geometry (precision-dependent, candidate generation only) ---------

    def _log_data(self, prec):
        """LLL-reduced log lattice of U^+ at the given precision."""
        R = RealField(prec)
        embs = self.K.real_embeddings(prec)
        rows = [[R(e(u)).log() for e in embs] for u in self.unit_basis]
        B = Matrix(R, rows)
        # LLL on the log lattice: reduce a rational approximation, apply the
        # transform to the exact generators.  Scaling factor is heuristic; a bad
        # reduction costs candidates, not correctness.
        scale = ZZ(2) ** (prec // 2)
        Bz = Matrix(ZZ, [[(x * scale).round() for x in row] for row in B.rows()])
        _, T = Bz.LLL(transformation=True)
        B = T.change_ring(R) * B
        units = []
        for row in T.rows():
            eps = self.K.one()
            for u, e in zip(self.unit_basis, row):
                eps *= u ** ZZ(e)
            units.append(eps)
        return R, embs, B, B * B.transpose(), units

    def _candidates(self, x, best_trace, prec):
        """
        Integer exponent vectors whose unit multiple *could* beat ``best_trace``.

        Uses the rigorous box:  if Tr(u) <= T then every log coordinate of u
        satisfies  log N(x) - (n-1) log T  <=  t_i  <=  log T.
        """
        R, embs, B, G, units = self._log_data(prec)
        n = self.n
        Lx = vector(R, [R(e(x)).log() for e in embs])
        logN = R(ZZ(self.K(x).norm()).abs()).log()

        M = R(best_trace).log()
        m = logN - (n - 1) * M
        centre = vector(R, [(m + M) / 2] * n)
        radius = R(n).sqrt() * (M - m) / 2

        tc = centre - Lx
        # y* solves G y = B tc; approximate solution is fine, radius is inflated.
        y = G.solve_right(B * tc)
        rho = radius * R(1.001) + R(2) ** (-prec // 4)

        exps = _enumerate_close(G, y, rho ** 2, self.rank, self.MAX_CANDIDATES)
        return units, exps

    # -- the main entry points --------------------------------------------

    def canonical(self, x, verify=True):
        """
        The canonical representative of the U^+-orbit of the totally positive ``x``.

        Raises ``ValueError`` if ``x`` is not totally positive.
        """
        x = self.K(x)
        if not x.is_totally_positive():
            raise ValueError(f"{x} is not totally positive")
        if self.rank == 0:                       # K = Q
            return x

        best = self._minimise(x, self.prec)
        if verify:
            again = self._minimise(x, 2 * self.prec)
            if again != best:
                # Low precision missed something; escalate until stable.
                prec = 4 * self.prec
                best = again
                while True:
                    nxt = self._minimise(x, prec)
                    if nxt == best:
                        break
                    best, prec = nxt, 2 * prec
                    if prec > 64 * self.prec:
                        raise RuntimeError("canonical representative did not stabilise")
        return best

    def _minimise(self, x, prec):
        """Minimal-trace element of the orbit, lex tie-break, at one precision."""
        best = x
        best_trace = ZZ(x.trace())
        best_coords = tuple(self.coordinates(x))

        while True:
            units, exps = self._candidates(x, best_trace, prec)
            improved = False
            for a in exps:
                eps = self.K.one()
                for u, e in zip(units, a):
                    if e:
                        eps *= u ** e
                y = eps * x
                t = ZZ(y.trace())
                if t > best_trace:
                    continue
                c = tuple(self.coordinates(y))
                if t < best_trace or c < best_coords:
                    if t < best_trace:
                        improved = True
                    best, best_trace, best_coords = y, t, c
            if not improved:
                return best

    def canonicalize(self, elements, verify=True):
        """Canonicalise a list, drop duplicate orbits, and sort deterministically."""
        seen = {}
        for x in elements:
            c = self.canonical(x, verify=verify)
            seen[tuple(self.coordinates(c))] = c
        return sorted(seen.values(), key=self.sort_key)

    def same_orbit(self, x, y):
        """Exact test: are ``x`` and ``y`` in the same U^+-orbit?"""
        x, y = self.K(x), self.K(y)
        if x.is_zero() or y.is_zero():
            return x.is_zero() and y.is_zero()
        if x.norm() != y.norm():
            return False
        q = x / y
        return q.is_integral() and (~q).is_integral() and q.is_totally_positive()

    def is_canonical(self, x):
        """Whether ``x`` is already the canonical representative of its orbit."""
        return self.K(x) == self.canonical(x)
