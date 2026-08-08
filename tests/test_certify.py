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


def report(failures, context):
    """
    Fail once, listing everything that went wrong.

    These tests loop over hundreds of elements.  Asserting inside the loop stops
    at the first bad one, which tells you nothing about whether it is an isolated
    case or the whole field -- and that distinction is usually the diagnosis.
    """
    if failures:
        head = "\n".join(f"    {f}" for f in failures[:15])
        more = f"\n    ... and {len(failures) - 15} more" if len(failures) > 15 else ""
        raise AssertionError(
            f"{context}: {len(failures)} mismatch(es)\n{head}{more}")

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


def _coordinate_box(ctx, trace_max):
    """
    Coordinate range that provably covers every totally positive element of
    trace at most ``trace_max``.

    A hardcoded box is not safe: it depends on how skewed the integral basis is.
    For Q(sqrt 10) with a +-8 box, every indecomposable of norm 10 -- including
    10 + 3*sqrt(10), which sits exactly at the Dress-Scharlau bound disc/4 --
    has first coordinate 10 and falls outside, so an entire norm class is
    invisible while the test still looks like it passes.

    If y is totally positive with Tr(y) <= T then 0 < sigma_j(y) < T for every
    j, and the coordinate vector is a = v * M^-1 with v the embedding vector and
    M the matrix of basis embeddings.  Hence |a_i| <= T * sum_j |M^-1[j][i]|.
    """
    n = ctx.degree
    embs = ctx.real_embeddings()
    M = sage.Matrix(sage.RDF, [[e(b) for e in embs] for b in ctx.basis])
    Minv = M.inverse()
    return [int(trace_max * sum(abs(Minv[j][i]) for j in range(n))) + 1
            for i in range(n)]


def _totally_positive_sample(ctx, trace_max=40):
    """Every totally positive element of trace at most ``trace_max``."""
    box = _coordinate_box(ctx, trace_max)
    out = []
    for coeffs in sage.cartesian_product([range(-b, b + 1) for b in box]):
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
    failures = []
    for y in sample:
        # Brute force over this sample is *exact*, because the sample is closed
        # downwards: if beta < y then beta is totally positive with
        # Tr(beta) < Tr(y) <= trace_max, so beta is in the sample too.  Hence
        # the comparison is an equality rather than a one-sided implication, and
        # a missed indecomposable is caught as well as a spurious one.
        brute = any((y - z).is_totally_positive() for z in lookup if z != y)
        got = is_indecomposable(y, ctx, use_norm_bound=False)
        if got != (not brute):
            failures.append(
                f"{y}  (norm {ZZ(y.norm())}, trace {ZZ(y.trace())}): "
                f"code says {'indecomposable' if got else 'decomposable'}, "
                f"exhaustive says {'decomposable' if brute else 'indecomposable'}")
    report(failures, f"Q(sqrt {D}), {len(sample)} elements of trace <= 30")


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
    sample = _totally_positive_sample(ctx, trace_max=60)
    found = {}
    for y in sample:
        if is_indecomposable(y, ctx):
            found.setdefault(ZZ(y.norm()), y)
    norms = set(found)
    detail = "  ".join(f"{n}:{found[n]}" for n in sorted(norms))
    assert norms == expected, (
        f"Q(sqrt {D}): got norms {sorted(norms)}, expected {sorted(expected)}\n"
        f"    sampled {len(sample)} elements of trace <= 60, "
        f"coordinate box {_coordinate_box(ctx, 60)}\n"
        f"    witnesses: {detail}")


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

def test_covering_radius_is_a_real_bound(cubic):
    """
    Babai's bound must actually bound: every log vector of a totally positive
    unit, reduced against the lattice, has to land within rho of the origin.

    This is the value that sits in the exponent of the search bound, so an
    error here silently truncates the enumeration rather than raising.
    """
    rho = cubic.log_unit_covering_radius
    assert rho > 0
    B = cubic.log_unit_lattice
    for row in B.rows():
        assert row.norm() > 0
    # the bound is at least half the longest reduced basis vector
    assert rho >= max(row.norm() for row in B.rows()) / 2 * 0.999


def test_trace_bound_exceeds_every_result(cubic):
    """Nothing found may exceed the bound used to find it."""
    from indecomposables.enumerate import indecomposables_exhaustive
    embs = cubic.real_embeddings()
    for y, _ in indecomposables_exhaustive(cubic):
        assert sum(embs[i](y) for i in range(cubic.degree)) <= cubic.trace_bound


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
def test_signature_argument_defaults_to_totally_positive(D):
    """Passing the all-plus signature explicitly must change nothing."""
    ctx = _ctx(D)
    plus = tuple([1] * ctx.degree)
    for y in _totally_positive_sample(ctx):
        assert (is_indecomposable(y, ctx, plus, use_norm_bound=False)
                == is_indecomposable(y, ctx, use_norm_bound=False))


@pytest.mark.parametrize("D", [5, 2, 13])
def test_s_indecomposability_in_every_signature_class(D):
    """
    The signed box test against a direct scan, in each class separately.

    Brute force here uses the *definition* -- is x a sum of two elements of the
    same class -- while the code uses the equivalent box condition, so this
    checks the equivalence as well as the implementation.
    """
    ctx = _ctx(D)
    embs = ctx.real_embeddings()
    box = _coordinate_box(ctx, 40)
    by_class = {}
    for coeffs in sage.cartesian_product([range(-b, b + 1) for b in box]):
        y = ctx.from_coordinates(coeffs)
        vals = [embs[i](y) for i in range(ctx.degree)]
        if any(v == 0 for v in vals) or sum(v ** 2 for v in vals) > 200:
            continue
        sig = tuple(1 if v > 0 else -1 for v in vals)
        by_class.setdefault(sig, []).append((y, vals))

    for sig, members in by_class.items():
        for y, vals in members:
            # Same downward-closure argument, in absolute value: anything
            # strictly inside the box of y is itself inside the T2 cutoff, so
            # this comparison is exact in both directions.
            others = {z for z, _ in members if z != y}
            brute = any((y - z) in others for z in others)
            assert is_indecomposable(y, ctx, sig) == (not brute), (
                f"{y} signature {sig}: exhaustive says "
                f"{'decomposable' if brute else 's-indecomposable'}")
