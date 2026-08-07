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
    "is_on_sail", "classify", "minimality_witness", "is_minimal",
]

#: Refuse to enumerate more than this many candidates before failing loudly.
MAX_CANDIDATES = 500000


# ---------------------------------------------------------------------------
# Lattice points in a box
# ---------------------------------------------------------------------------

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
    y = G.solve_right(M * centre)
    rho = R(n).sqrt() / 2 * R(1.001) + R(2) ** (-prec // 4)

    return _enumerate_close(G, y, rho ** 2, len(basis), cap)


# ---------------------------------------------------------------------------
# Decomposability
# ---------------------------------------------------------------------------

def minimality_witness(x, ctx, signature=None, prec=None):
    """
    A ``beta`` in the same signature class with ``|sigma_i(beta)| < |sigma_i(x)|``
    for every ``i``, or ``None``.

    For the totally positive class this is exactly a decomposition: the
    condition reads ``0 < beta < x``, so ``x = beta + (x - beta)`` with both
    summands totally positive.  For other classes ``S_sigma`` is not closed
    under addition, so no decomposition is implied and the right word is
    *minimal* rather than *indecomposable* -- but the geometry, and hence the
    computation, is identical.
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


def is_minimal(x, ctx, signature=None, verify=True):
    """Whether ``x`` is minimal in its signature class."""
    if minimality_witness(x, ctx, signature, ctx.prec) is not None:
        return False
    if verify and minimality_witness(x, ctx, signature, 2 * ctx.prec) is not None:
        return False
    return True


def decomposition_witness(x, ctx, prec=None):
    """
    A ``beta`` with ``0 < beta < x``, or ``None`` if the search found none.

    A returned witness is exact and needs no further checking: it is verified
    with :meth:`is_totally_positive` before being handed back.  A ``None`` is
    only as complete as the working precision, which is why
    :func:`is_indecomposable`, not this function, is the public test.
    """
    K = ctx.K
    x = K(x)
    if not x.is_totally_positive():
        raise ValueError(f"{x} is not totally positive")
    prec = prec or ctx.prec

    embs = ctx.real_embeddings(prec)
    uppers = [embs[i](x) for i in range(ctx.degree)]
    for a in _box_points(ctx.basis, uppers, ctx, prec):
        beta = ctx.from_coordinates(a)
        if beta.is_zero() or beta == x:
            continue
        if beta.is_totally_positive() and (x - beta).is_totally_positive():
            return beta
    return None


def is_indecomposable(x, ctx, verify=True, use_norm_bound=True):
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
    if not x.is_totally_positive():
        raise ValueError(f"{x} is not totally positive")

    if use_norm_bound and abs(ZZ(x.norm())) > abs(ctx.discriminant):
        return False

    if codifferent_certificate(x, ctx) is not None:
        return True

    if decomposition_witness(x, ctx, ctx.prec) is not None:
        return False
    if verify and decomposition_witness(x, ctx, 2 * ctx.prec) is not None:
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

    beta = decomposition_witness(x, ctx, ctx.prec)
    if beta is None and verify:
        beta = decomposition_witness(x, ctx, 2 * ctx.prec)
    if beta is not None:
        return False, False, "box_witness", beta
    return True, False, "box_exhausted", None
