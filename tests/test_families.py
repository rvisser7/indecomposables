"""
Cross-validation: every specialised algorithm against exhaustive search.

This is the most valuable test in the project, and the one unit tests
structurally cannot replace.  A hand-written classification and a hand-written
test can encode the *same* misunderstanding of a theorem; brute force cannot,
because it computes the definition directly.

So: for every field where a specialised algorithm applies, run it and run the
exhaustive method, canonicalise both, and demand they agree exactly.  The
comparison is on canonical representatives, which is why pinning the canonical
form mattered -- without it the two could return equivalent-but-different lists
and the test would fail for no mathematical reason.

Marked ``verification`` because it is slow.  Run it before cutting a data
release, not on every push::

    sage -python -m pytest tests/test_families.py -m verification
"""

import pytest

sage = pytest.importorskip("sage.all", reason="Sage not available")
NumberField = sage.NumberField
PolynomialRing = sage.PolynomialRing
QQ, ZZ = sage.QQ, sage.ZZ

pytestmark = [pytest.mark.sage, pytest.mark.comparison]

from indecomposables.context import FieldContext                    # noqa: E402
from indecomposables.enumerate import indecomposables_exhaustive    # noqa: E402
from indecomposables.record import applicable_algorithms            # noqa: E402

R = PolynomialRing(QQ, "x")
x = R.gen()


def _canonical_set(ctx, pairs):
    """Compare on canonical coordinate vectors, not on element objects."""
    return {tuple(ctx.coordinates(ctx.normalizer.canonical(y))) for y, _ in pairs}


def _sail_set(ctx, pairs):
    return {tuple(ctx.coordinates(ctx.normalizer.canonical(y)))
            for y, on_sail in pairs if on_sail}


# ---------------------------------------------------------------------------
# Field corpora
# ---------------------------------------------------------------------------

def real_quadratic(dmax=60):
    """Squarefree D > 1, as x^2 - x - (D-1)/4 or x^2 - D."""
    out = []
    for D in range(2, dmax):
        if any(D % (p * p) == 0 for p in range(2, int(D ** 0.5) + 1)):
            continue
        f = x ** 2 - x - (D - 1) // 4 if D % 4 == 1 else x ** 2 - D
        out.append((f"Q(sqrt {D})", f))
    return out


def simplest_cubics(nmax=25):
    """x^3 - n x^2 - (n+3) x - 1, the simplest cubic family."""
    return [(f"simplest cubic n={n}", x ** 3 - n * x ** 2 - (n + 3) * x - 1)
            for n in range(-1, nmax)]


# ---------------------------------------------------------------------------
# The cross-check
# ---------------------------------------------------------------------------

def _cross_check(name, f):
    K = NumberField(f, "a")
    if not K.is_totally_real():
        pytest.skip(f"{name} is not totally real")
    ctx = FieldContext(K, label=name)

    algorithms = applicable_algorithms(ctx)
    specialised = [(n, fn) for n, fn in algorithms
                   if n not in ("brute_force", "sail_walk")]
    if not specialised:
        pytest.skip(f"no specialised algorithm applies to {name}")

    reference = indecomposables_exhaustive(ctx)
    want = _canonical_set(ctx, reference)

    for algo_name, fn in specialised:
        got = _canonical_set(ctx, fn(ctx))
        missing, extra = want - got, got - want
        assert not missing, (
            f"{name}: {algo_name} MISSED {len(missing)} indecomposable(s) "
            f"that brute force found, e.g. {sorted(missing)[:3]}")
        assert not extra, (
            f"{name}: {algo_name} returned {len(extra)} element(s) that are NOT "
            f"indecomposable, e.g. {sorted(extra)[:3]}")


@pytest.mark.verification
@pytest.mark.dress_scharlau
@pytest.mark.parametrize("name,f", real_quadratic())
def test_real_quadratic_matches_brute_force(name, f):
    _cross_check(name, f)


@pytest.mark.verification
@pytest.mark.kala_tinkova
@pytest.mark.parametrize("name,f", simplest_cubics())
def test_simplest_cubic_matches_brute_force(name, f):
    _cross_check(name, f)


# ---------------------------------------------------------------------------
# Sail vs indecomposable: the mathematically interesting comparison
# ---------------------------------------------------------------------------

@pytest.mark.verification
@pytest.mark.sail
@pytest.mark.parametrize("name,f", real_quadratic(40))
def test_degree_two_sail_is_everything(name, f):
    """Degree 2 is the case where sail and indecomposable coincide."""
    K = NumberField(f, "a")
    ctx = FieldContext(K, label=name)
    pairs = indecomposables_exhaustive(ctx)
    assert _sail_set(ctx, pairs) == _canonical_set(ctx, pairs)


@pytest.mark.verification
@pytest.mark.sail
@pytest.mark.parametrize("name,f", simplest_cubics(15))
def test_sail_is_contained_in_indecomposables(name, f):
    """
    Degree >= 3: the sail is a subset, possibly proper.  The size of the gap is
    the novel measurement the database is for, so this records it rather than
    asserting it is empty.
    """
    K = NumberField(f, "a")
    if not K.is_totally_real():
        pytest.skip("not totally real")
    ctx = FieldContext(K, label=name)
    pairs = indecomposables_exhaustive(ctx)
    all_set, sail = _canonical_set(ctx, pairs), _sail_set(ctx, pairs)
    assert sail <= all_set, f"{name}: a sail point is not indecomposable"
    print(f"{name}: {len(all_set)} indecomposable, {len(sail)} on the sail, "
          f"{len(all_set) - len(sail)} off it")


# ---------------------------------------------------------------------------
# Contract conformance -- fast, so it runs in the normal tier
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,f", real_quadratic(20) + simplest_cubics(6))
def test_family_returns_the_declared_shape(name, f):
    """
    Every registered algorithm must return ``[(element, on_sail), ...]`` with
    totally positive elements.  Cheap, and it catches an interface drift before
    the slow cross-check ever runs.
    """
    K = NumberField(f, "a")
    if not K.is_totally_real():
        pytest.skip("not totally real")
    ctx = FieldContext(K, label=name)
    for algo_name, fn in applicable_algorithms(ctx):
        if algo_name == "brute_force":
            continue
        out = fn(ctx)
        assert isinstance(out, (list, tuple)), f"{algo_name} did not return a list"
        for item in out:
            assert len(item) == 2, f"{algo_name} items must be (element, on_sail)"
            y, on_sail = item
            assert ctx.K(y).is_totally_positive(), \
                f"{algo_name} returned a non totally positive element"
            assert on_sail in (True, False, 0, 1)
