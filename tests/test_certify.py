"""
Decomposability tests.

The interesting assertions here are relations that must hold for any input,
plus a handful of values pinned from Dress-Scharlau in the real quadratic case
where the classification is known.
"""

import pytest

sage = pytest.importorskip("sage.all", reason="Sage not available")
NumberField = sage.NumberField
PolynomialRing = sage.PolynomialRing
QQ, ZZ = sage.QQ, sage.ZZ

pytestmark = pytest.mark.sage

from indecomposables.context import FieldContext          # noqa: E402
from indecomposables.certify import (                     # noqa: E402
    classify, codifferent_certificate, decomposition_witness,
    is_indecomposable, is_on_sail,
)

R = PolynomialRing(QQ, "x")
x = R.gen()

# Norms of the indecomposables, up to totally positive units, for Q(sqrt D).
# Independently reproduced by exhaustive search; consistent with Dress-Scharlau
# and with the classical bound N(alpha) <= disc/4 in degree 2.
QUADRATIC_NORMS = {
    5: {1},
    2: {1, 2},
    3: {1},
    13: {1, 3},
    6: {1, 3},
    7: {1, 2},
    10: {1, 6, 9, 10},
    21: {1},
    33: {1, 3, 4},
}


def _ctx(D):
    f = x ** 2 - x - (D - 1) // 4 if D % 4 == 1 else x ** 2 - D
    return FieldContext(NumberField(f, "a"))


@pytest.fixture(scope="session")
def cubic():
    return FieldContext(NumberField(x ** 3 - x ** 2 - 2 * x + 1, "a"))


def _totally_positive_sample(ctx, trace_max=40):
    out = []
    for coeffs in sage.cartesian_product([range(-8, 9)] * ctx.degree):
        y = ctx.from_coordinates(coeffs)
        if y.is_totally_positive() and ZZ(y.trace()) <= trace_max:
            out.append(y)
    return out


# ---------------------------------------------------------------------------
# The definition itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", [5, 2, 3, 13, 6])
def test_witness_is_a_real_decomposition(D):
    ctx = _ctx(D)
    for y in _totally_positive_sample(ctx):
        beta = decomposition_witness(y, ctx)
        if beta is not None:
            assert beta.is_totally_positive()
            assert (y - beta).is_totally_positive()
            assert beta + (y - beta) == y


@pytest.mark.parametrize("D", [5, 2, 13, 7])
def test_agrees_with_exhaustive_search(D):
    """The box method against a direct scan of the coefficient grid."""
    ctx = _ctx(D)
    sample = _totally_positive_sample(ctx, trace_max=30)
    lookup = set(sample)
    for y in sample:
        brute = any((y - z).is_totally_positive() for z in lookup if z != y)
        assert is_indecomposable(y, ctx) != brute or not brute


def test_one_is_always_indecomposable(cubic):
    assert is_indecomposable(cubic.K.one(), cubic)


def test_rejects_non_totally_positive(cubic):
    with pytest.raises(ValueError):
        is_indecomposable(cubic.K(-1), cubic)


# ---------------------------------------------------------------------------
# Certificates
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", [5, 2, 13, 6, 10])
def test_certificate_implies_indecomposable(D):
    """A lambda with Tr(lambda x) = 1 settles it, with no completeness claim."""
    ctx = _ctx(D)
    for y in _totally_positive_sample(ctx):
        lam = codifferent_certificate(y, ctx)
        if lam is not None:
            assert lam.is_totally_positive()
            assert QQ((lam * y).trace()) == 1
            assert is_indecomposable(y, ctx)


@pytest.mark.parametrize("D", [5, 2, 13, 6, 10, 33])
def test_degree_two_sail_equals_indecomposable(D):
    """In degree 2 every indecomposable lies on the sail; the converse is general."""
    ctx = _ctx(D)
    for y in _totally_positive_sample(ctx):
        assert is_on_sail(y, ctx) == is_indecomposable(y, ctx)


# ---------------------------------------------------------------------------
# Known values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D,expected", sorted(QUADRATIC_NORMS.items()))
def test_quadratic_indecomposable_norms(D, expected):
    ctx = _ctx(D)
    norms = {ZZ(y.norm()) for y in _totally_positive_sample(ctx, trace_max=60)
             if is_indecomposable(y, ctx)}
    assert norms == expected


@pytest.mark.parametrize("D", sorted(QUADRATIC_NORMS))
def test_dress_scharlau_norm_bound(D):
    """In degree 2 indecomposables satisfy N(alpha) <= disc/4, sharper than |disc|."""
    ctx = _ctx(D)
    disc = abs(ZZ(ctx.discriminant))
    for y in _totally_positive_sample(ctx, trace_max=60):
        if is_indecomposable(y, ctx):
            assert 4 * ZZ(y.norm()) <= disc


# ---------------------------------------------------------------------------
# classify() bookkeeping
# ---------------------------------------------------------------------------

def test_classify_reasons_are_consistent(cubic):
    seen = set()
    for y in _totally_positive_sample(cubic, trace_max=25):
        indec, on_sail, reason, witness = classify(y, cubic)
        seen.add(reason)
        assert reason in {"norm_bound", "codifferent", "box_witness", "box_exhausted"}
        if reason == "codifferent":
            assert indec and on_sail and witness is not None
        if reason == "box_witness":
            assert not indec and witness.is_totally_positive()
        if reason == "norm_bound":
            assert not indec
            assert abs(ZZ(y.norm())) > abs(ZZ(cubic.discriminant))
        assert on_sail <= indec          # on the sail implies indecomposable
    assert "codifferent" in seen


@pytest.mark.slow
def test_precision_does_not_change_the_answer(cubic):
    lo = FieldContext(cubic.K, basis=cubic.basis, prec=32)
    hi = FieldContext(cubic.K, basis=cubic.basis, prec=256)
    for y in _totally_positive_sample(cubic, trace_max=25):
        assert is_indecomposable(y, lo) == is_indecomposable(y, hi)


# ---------------------------------------------------------------------------
# Signature classes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", [5, 2, 13, 6])
def test_minimality_generalises_decomposability(D):
    """
    On the totally positive class, minimal and indecomposable must coincide --
    the box condition |sigma_i(beta)| < |sigma_i(alpha)| reads 0 < beta < alpha
    there.
    """
    from indecomposables.certify import is_minimal
    ctx = _ctx(D)
    for y in _totally_positive_sample(ctx):
        assert is_minimal(y, ctx) == is_indecomposable(y, ctx, use_norm_bound=False)


@pytest.mark.parametrize("D", [5, 2, 13])
def test_minimality_in_every_signature_class(D):
    """The signed box test against a direct scan, in each class separately."""
    from indecomposables.certify import is_minimal
    ctx = _ctx(D)
    embs = ctx.real_embeddings()
    by_class = {}
    for coeffs in sage.cartesian_product([range(-10, 11)] * ctx.degree):
        y = ctx.from_coordinates(coeffs)
        vals = [embs[i](y) for i in range(ctx.degree)]
        if any(v == 0 for v in vals) or sum(v ** 2 for v in vals) > 200:
            continue
        sig = tuple(1 if v > 0 else -1 for v in vals)
        by_class.setdefault(sig, []).append((y, vals))

    for sig, members in by_class.items():
        for y, vals in members:
            brute = any(all(abs(w[i]) < abs(vals[i]) for i in range(ctx.degree))
                        for z, w in members if z != y)
            assert is_minimal(y, ctx, sig) != brute or not brute
