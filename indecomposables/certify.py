"""
Deciding indecomposability, and certifying it.

The definition is directly checkable: a totally positive ``alpha`` is
**decomposable** exactly when some ``beta`` in the order satisfies
``0 < beta < alpha`` in the totally positive order, since then
``alpha = beta + (alpha - beta)`` with both summands totally positive.  Under
the Minkowski embedding that asks for a lattice point in the open box

    prod_i (0, sigma_i(alpha)),

whose volume is ``N(alpha)`` against a covolume of ``sqrt(|disc|)``.  With the
Kala-Yatsyna bound ``N(alpha) <= |disc|`` for indecomposables, the expected
number of points is at most ``sqrt(|disc|)`` -- so the test is cheap, and it
early-exits on the first witness.

No mixed-integer program is involved.  A solver would be settling a yes/no
mathematical question through floating point; here floating point only ever
*proposes* candidates, and every accepted answer is checked exactly with
``is_totally_positive``.  A decomposition returned by this module is therefore
certainly a decomposition.  The other direction -- claiming indecomposability --
rests on the candidate set being complete, so :func:`is_indecomposable` repeats
the search at higher precision to confirm, as :mod:`indecomposables.normalize`
does.

Two cheap sufficient tests come first, in cost order:

* ``N(alpha) > |disc|`` implies decomposable (Kala-Yatsyna);
* a totally positive ``lambda`` in the codifferent with ``Tr(lambda alpha) = 1``
  implies indecomposable, since ``Tr(lambda beta) >= 1`` for every totally
  positive ``beta`` in the order, so a decomposition would give ``1 >= 2``.
  Such a ``lambda`` is a supporting hyperplane of the sail at ``alpha``, so the
  same computation decides sail membership.

Both searches are lattice points in a box -- the second in
``prod_i (0, 1/sigma_i(alpha))`` for the codifferent, because ``lambda`` and
``alpha`` are both totally positive with ``Tr(lambda alpha) = 1``, forcing every
``sigma_i(lambda)`` into ``(0, 1/sigma_i(alpha))``.  So one enumerator serves
both.
"""

from __future__ import annotations

from sage.all import Matrix, QQ, RealField, ZZ, vector

from .normalize import _enumerate_close

__all__ = [
    "decomposition_witness", "is_indecomposable", "codifferent_certificate",
    "codifferent_trace", "is_on_sail", "classify",
]

#: Refuse to enumerate more than this many candidates before failing loudly.
MAX_CANDIDATES = 500000


# ---------------------------------------------------------------------------
# Lattice points in a box
# ---------------------------------------------------------------------------

def _solve_spd(G, rhs):
    """
    Solve ``G y = rhs`` for symmetric positive definite ``G``, by Cholesky.

    Written out rather than delegated because Sage's linear algebra over
    ``RealField(prec)`` is patchy -- ``gram_schmidt`` is simply unimplemented
    there -- and dropping to ``RDF`` would cap the working precision at 53 bits.
    """
    n = G.nrows()
    R = G.base_ring()
    L = [[R(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            acc = G[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if acc <= 0:
                    raise ValueError(
                        "Gram matrix is not positive definite at this precision")
                L[i][i] = acc.sqrt()
            else:
                L[i][j] = acc / L[j][j]
    # forward then back substitution
    z = [R(0)] * n
    for i in range(n):
        z[i] = (rhs[i] - sum(L[i][k] * z[k] for k in range(i))) / L[i][i]
    y = [R(0)] * n
    for i in reversed(range(n)):
        y[i] = (z[i] - sum(L[k][i] * y[k] for k in range(i + 1, n))) / L[i][i]
    return vector(R, y)

def _box_points(basis, uppers, ctx, prec, signs=None, cap=MAX_CANDIDATES):
    """
    Integer coefficient vectors ``a`` with ``0 < s_i (a . basis)_i < uppers_i``.

    ``signs`` is a vector of +-1 selecting an orthant, defaulting to the totally
    positive one.  Flipping a sign reflects that coordinate, which is why the
    same routine serves every signature class.

    ``basis`` is a Z-basis of the lattice as field elements and ``uppers`` are
    positive reals.  Rescaling coordinate ``i`` by ``1 / uppers_i`` turns the box
    into the unit cube, which sits inside the ball of radius ``sqrt(n)/2`` about
    its centre -- so the search reduces to the close-vector enumeration already
    used for canonical representatives.

    Returns a *superset*; the caller filters exactly.
    """
    R = RealField(prec)
    n = ctx.degree
    embs = ctx.real_embeddings(prec)

    if signs is None:
        signs = [1] * n
    M = Matrix(R, [[R(s) * R(e(b)) / R(u)
                    for e, u, s in zip(embs, uppers, signs)] for b in basis])
    G = M * M.transpose()
    centre = vector(R, [R(1) / 2] * n)
    # G is symmetric positive definite, so a Cholesky solve works and avoids
    # relying on Sage's generic solver over an inexact ring -- the same class of
    # gap as `gram_schmidt`, which is unimplemented for RealField.
    y = _solve_spd(G, M * centre)
    rho = R(n).sqrt() / 2 * R(1.001) + R(2) ** (-prec // 4)

    return _enumerate_close(G, y, rho ** 2, len(basis), cap)


# ---------------------------------------------------------------------------
# Decomposability
# ---------------------------------------------------------------------------

def decomposition_witness(x, ctx, signature=None, prec=None):
    """
    A ``beta`` in the same signature class with ``x = beta + gamma``, or ``None``.

    ``x`` in the signature class ``S_s`` is **s-indecomposable** when it is not a
    sum of two elements of ``S_s``.  Each ``S_s`` is the lattice intersected with
    an open orthant, hence closed under addition, so this is the same notion in
    every class; for ``s`` totally positive it is ordinary indecomposability.

    The search uses the equivalent geometric form -- some ``beta`` in ``S_s``
    with ``|sigma_i(beta)| < |sigma_i(x)|`` for every ``i`` -- because that is a
    box.  The two agree: a decomposition gives
    ``|sigma_i(x)| = |sigma_i(beta)| + |sigma_i(gamma)| > |sigma_i(beta)|``, and
    conversely ``gamma = x - beta`` keeps every sign ``s_i`` and is nonzero.

    A returned witness is exact and needs no further checking: both it and
    ``x - beta`` are verified to lie in ``S_s`` before it is handed back.
    """
    K = ctx.K
    x = K(x)
    prec = prec or ctx.prec
    embs = ctx.real_embeddings(prec)
    vals = [embs[i](x) for i in range(ctx.degree)]
    if signature is None:
        signature = [1 if v > 0 else -1 for v in vals]
    if any(s * v <= 0 for s, v in zip(signature, vals)):
        raise ValueError(f"{x} does not have signature {tuple(signature)}")

    uppers = [abs(v) for v in vals]
    for a in _box_points(ctx.basis, uppers, ctx, prec, signs=signature):
        beta = ctx.from_coordinates(a)
        if beta.is_zero() or beta == x:
            continue
        bvals = [embs[i](beta) for i in range(ctx.degree)]
        if all(s * b > 0 and abs(b) < abs(v)
               for s, b, v in zip(signature, bvals, vals)):
            return beta
    return None


def is_indecomposable(x, ctx, signature=None, verify=True, use_norm_bound=True):
    """
    Whether ``x`` is indecomposable.

    With ``verify`` the search is repeated at doubled precision before ``True``
    is returned, since that answer rests on the candidate set being complete.
    ``False`` always comes with an exact witness and needs no confirmation.

    EXAMPLES::

        sage: from indecomposables import FieldContext           # doctest: +SKIP
        sage: from indecomposables.certify import is_indecomposable
        sage: ctx = FieldContext.from_coeffs([-1, -1, 1])         # doctest: +SKIP
        sage: is_indecomposable(ctx.K.one(), ctx)                 # doctest: +SKIP
        True
    """
    K = ctx.K
    x = K(x)
    totally_positive = signature is None or all(s > 0 for s in signature)

    if totally_positive:
        # Two cheap sufficient tests, available only in the totally positive
        # class: the Kala-Yatsyna norm bound, and a codifferent certificate.
        if use_norm_bound and abs(ZZ(x.norm())) > abs(ctx.discriminant):
            return False
        if codifferent_certificate(x, ctx) is not None:
            return True

    if decomposition_witness(x, ctx, signature, ctx.prec) is not None:
        return False
    if verify and decomposition_witness(x, ctx, signature, 2 * ctx.prec) is not None:
        return False
    return True


# ---------------------------------------------------------------------------
# Certificates and the sail
# ---------------------------------------------------------------------------

def _codifferent_basis(ctx):
    if getattr(ctx, "_codiff_basis", None) is None:
        ctx._codiff_basis = list((ctx.K.different() ** -1).basis())
    return ctx._codiff_basis


def codifferent_certificate(x, ctx, prec=None):
    """
    A totally positive ``lambda`` in the codifferent with ``Tr(lambda x) = 1``.

    Such a ``lambda`` proves ``x`` indecomposable outright, with no completeness
    assumption, and is the normal vector of a supporting hyperplane of the sail
    at ``x``.

    ``None`` does not mean ``x`` is decomposable when the degree is at least 3:
    it means ``x`` is either decomposable or an indecomposable lying off the
    sail, and the box test has to decide which.
    """
    K = ctx.K
    x = K(x)
    prec = prec or ctx.prec
    embs = ctx.real_embeddings(prec)
    uppers = [1 / embs[i](x) for i in range(ctx.degree)]
    basis = _codifferent_basis(ctx)

    for a in _box_points(basis, uppers, ctx, prec):
        lam = sum(ZZ(c) * b for c, b in zip(a, basis))
        if lam.is_zero() or not lam.is_totally_positive():
            continue
        if QQ((lam * x).trace()) == 1:
            return lam
    return None


def codifferent_trace(x, ctx, prec=None):
    """
    ``min { Tr(lambda x) : lambda > 0 in the codifferent }``, a positive integer.

    Equals 1 exactly when ``x`` lies on the sail, so this is a graded refinement
    of :func:`is_on_sail`: it says not merely that ``x`` is off the sail but how
    far.  Known to be 1 throughout for real quadratic fields and at most 2 for
    simplest cubics; Kala-Tinkova show it can be arbitrarily large in other
    families, and its distribution is not known -- which is the reason to record
    it across a large database.

    Costs nothing beyond :func:`codifferent_certificate`, which enumerates the
    same box and stops at 1 instead of taking the minimum.
    """
    K = ctx.K
    x = K(x)
    prec = prec or ctx.prec
    embs = ctx.real_embeddings(prec)
    n = ctx.degree
    basis = _codifferent_basis(ctx)

    # Tr(lambda x) = t with lambda, x both totally positive forces
    # 0 < sigma_i(lambda) < t / sigma_i(x).  The bound is unknown in advance, so
    # widen the box until something is found; in practice the first try succeeds.
    best = None
    for t in (1, 2, 4, 8, 16, 32):
        uppers = [t / embs[i](x) for i in range(n)]
        for a in _box_points(basis, uppers, ctx, prec):
            lam = sum(ZZ(c) * b for c, b in zip(a, basis))
            if lam.is_zero() or not lam.is_totally_positive():
                continue
            tr = QQ((lam * x).trace())
            if tr > 0 and tr in ZZ and (best is None or tr < best):
                best = ZZ(tr)
        if best is not None:
            return best
    return None


def is_on_sail(x, ctx):
    """
    Whether ``x`` lies on the sail.

    In degree 2 this coincides with indecomposability.  In degree at least 3 the
    sail's lattice points are all indecomposable but the converse fails, so this
    is strictly stronger, and the gap between the two is itself worth recording.
    """
    return codifferent_certificate(x, ctx) is not None


def classify(x, ctx, verify=True):
    """
    Decide ``x``, and record which test decided it.

    Returns ``(indecomposable, on_sail, reason, witness)`` with ``reason`` one of
    ``norm_bound``, ``codifferent``, ``box_witness``, ``box_exhausted``.  Rows
    settled by different routes carry different assumptions -- only
    ``box_exhausted`` depends on completeness of a numerical search -- so
    keeping the reason makes the database auditable after the fact.
    """
    K = ctx.K
    x = K(x)
    if not x.is_totally_positive():
        raise ValueError(f"{x} is not totally positive")

    if abs(ZZ(x.norm())) > abs(ctx.discriminant):
        return False, False, "norm_bound", None

    lam = codifferent_certificate(x, ctx)
    if lam is not None:
        return True, True, "codifferent", lam

    beta = decomposition_witness(x, ctx, None, ctx.prec)
    if beta is None and verify:
        beta = decomposition_witness(x, ctx, None, 2 * ctx.prec)
    if beta is not None:
        return False, False, "box_witness", beta
    return True, False, "box_exhausted", None
