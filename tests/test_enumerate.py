"""
Enumeration: bounds, candidate generation, and the full per-field computation.
"""

import pytest

sage = pytest.importorskip("sage.all", reason="Sage not available")
NumberField = sage.NumberField
PolynomialRing = sage.PolynomialRing
QQ, ZZ = sage.QQ, sage.ZZ

pytestmark = pytest.mark.sage

from indecomposables.context import FieldContext                    # noqa: E402
from indecomposables.certify import is_indecomposable               # noqa: E402
from indecomposables.enumerate import (                             # noqa: E402
    all_signature_classes, brunotte_norm_bound, candidates,
    candidates_by_norm, indecomposables_exhaustive, indecomposables_in_class,
    norm_bound_for, t2_bound, trace_form, unit_representatives,
)

R = PolynomialRing(QQ, "x")
x = R.gen()

FIELDS = {
    "2.2.5.1": x ** 2 - x - 1,
    "2.2.8.1": x ** 2 - 2,
    "2.2.40.1": x ** 2 - 10,          # indecomposable norms {1, 6, 9, 10}
    "3.3.49.1": x ** 3 - x ** 2 - 2 * x + 1,
    "3.3.81.1": x ** 3 - 3 * x - 1,
    "4.4.725.1": x ** 4 - x ** 3 - 3 * x ** 2 + x + 1,
}


@pytest.fixture(scope="session")
def ctxs():
    return {lab: FieldContext(NumberField(f, "a"), label=lab)
            for lab, f in FIELDS.items()}


@pytest.fixture(params=sorted(FIELDS))
def ctx(request, ctxs):
    return ctxs[request.param]


def test_trace_form_is_positive_definite(ctx):
    G = trace_form(ctx)
    assert G.is_symmetric()
    assert all(G[:k, :k].determinant() > 0 for k in range(1, ctx.degree + 1))


def test_t2_bound_covers_every_result(ctx):
    """
    Nothing found may exceed the bound that was used to find it -- otherwise the
    search truncated and the answer is incomplete.
    """
    B = t2_bound(ctx)
    embs = ctx.real_embeddings()
    for y, _ in indecomposables_exhaustive(ctx):
        t2 = sum(embs[i](y) ** 2 for i in range(ctx.degree))
        assert t2 <= B


def test_candidates_respect_the_bound(ctx):
    """
    candidates() returns a superset: the bound is rounded up to an integer for
    PARI, so a vector may exceed it by less than 1.  Rounding down instead would
    risk dropping a vector sitting exactly on the bound, which is the wrong
    trade -- the extras are removed by the norm filter anyway.
    """
    from indecomposables.enumerate import CandidateExplosion
    B = t2_bound(ctx)
    limit = int(B) + 1
    embs = ctx.real_embeddings()
    n = 0
    try:
        for v in candidates(ctx, B):
            y = ctx.from_coordinates(v)
            assert sum(embs[i](y) ** 2 for i in range(ctx.degree)) <= limit
            n += 1
            if n > 500:
                break
    except CandidateExplosion:
        pytest.skip("lattice enumeration too large for this field; "
                    "the ideal method is the one used in anger")
    assert n > 0


def test_everything_returned_is_indecomposable(ctx):
    for y, _ in indecomposables_exhaustive(ctx):
        assert y.is_totally_positive()
        assert is_indecomposable(y, ctx)


def test_one_is_present(ctx):
    coords = {tuple(ctx.coordinates(y)) for y, _ in indecomposables_exhaustive(ctx)}
    assert tuple(ctx.coordinates(ctx.K.one())) in coords


def test_results_are_canonical_and_deduplicated(ctx):
    out = [y for y, _ in indecomposables_exhaustive(ctx)]
    coords = [tuple(ctx.coordinates(y)) for y in out]
    assert len(coords) == len(set(coords))
    for y in out:
        assert ctx.normalizer.is_canonical(y)


def test_norms_respect_kala_yatsyna(ctx):
    disc = abs(ZZ(ctx.discriminant))
    for y, _ in indecomposables_exhaustive(ctx):
        assert abs(ZZ(y.norm())) <= disc


def test_sail_is_a_subset(ctx):
    """Sail points are always indecomposable; the converse fails for degree >= 3."""
    out = indecomposables_exhaustive(ctx)
    assert all(is_indecomposable(y, ctx) for y, on_sail in out if on_sail)
    if ctx.degree == 2:
        assert all(on_sail for _, on_sail in out)


def test_signature_class_zero_is_totally_positive(ctx):
    masks, sigs, results = all_signature_classes(ctx)
    assert masks[0] == 0
    assert sigs[0] == tuple([1] * ctx.degree)
    tp = {tuple(ctx.coordinates(y)) for y, _, _ in results[0]}
    assert tp == {tuple(ctx.coordinates(y)) for y, _ in indecomposables_exhaustive(ctx)}


def test_signature_class_count_matches_unit_rank(ctx):
    masks, sigs, results = all_signature_classes(ctx)
    assert len(masks) == 2 ** (ctx.degree - ctx.unit_signature_rank)
    assert len(results) == len(sigs)


def test_elements_lie_in_their_declared_class(ctx):
    masks, sigs, results = all_signature_classes(ctx)
    embs = ctx.real_embeddings()
    for sig, res in zip(sigs, results):
        for y, _, _ in res:
            got = tuple(1 if embs[i](y) > 0 else -1 for i in range(ctx.degree))
            assert got == tuple(sig)


@pytest.mark.parametrize("label,expected", [
    ("2.2.5.1", {1}),
    ("2.2.8.1", {1, 2}),
    ("2.2.40.1", {1, 6, 9, 10}),
])
def test_known_quadratic_norms(label, expected, ctxs):
    ctx = ctxs[label]
    norms = {abs(ZZ(y.norm())) for y, _ in indecomposables_exhaustive(ctx)}
    assert norms == expected


@pytest.mark.slow
def test_norm_bound_is_not_binding(ctx):
    """Doubling the norm bound must not find anything new."""
    base = {tuple(ctx.coordinates(y)) for y, _, _ in indecomposables_in_class(ctx)}
    wider = {tuple(ctx.coordinates(y))
             for y, _, _ in indecomposables_in_class(
                 ctx, bound=2 * abs(ZZ(ctx.discriminant)))}
    assert base == wider


# ---------------------------------------------------------------------------
# The two candidate generators
# ---------------------------------------------------------------------------

def test_ideal_and_lattice_methods_agree(ctx):
    """
    The most valuable test here.

    The two generators share nothing: one enumerates principal ideals of bounded
    norm, the other lattice points of bounded T2 via qfminim.  They must find
    exactly the same s-indecomposables, so this checks both the ideal sweep and
    the T2 bound at once -- a bug in either shows up as a difference.
    """
    from indecomposables.enumerate import CandidateExplosion
    by_ideal = {tuple(ctx.coordinates(y))
                for y, _, _ in indecomposables_in_class(ctx, method="ideals")}
    try:
        by_lattice = {tuple(ctx.coordinates(y))
                      for y, _, _ in indecomposables_in_class(ctx, method="lattice")}
    except CandidateExplosion as exc:
        # Expected from degree 4 up: the T2 bound carries a factor e^(n rho)
        # that ideal enumeration does not pay.  Not a failure -- it is the
        # reason the ideal method is the default.
        pytest.skip(f"lattice method not feasible here: {str(exc)[:120]}")
    missing, extra = by_lattice - by_ideal, by_ideal - by_lattice
    assert not missing, f"the ideal sweep missed {sorted(missing)[:3]}"
    assert not extra, f"the ideal sweep found extra {sorted(extra)[:3]}"


def test_unit_representatives_cover_every_signature(ctx):
    """
    One sweep over ideals has to reach every signature class, which it does only
    if the representatives cover U/U^2.
    """
    embs = ctx.real_embeddings()
    sigs = set()
    for u in unit_representatives(ctx):
        sigs.add(tuple(1 if embs[i](u) > 0 else -1 for i in range(ctx.degree)))
    assert len(sigs) == 2 ** ctx.unit_signature_rank


def test_norm_bound_is_the_better_of_the_two(ctx):
    """Both bounds are proved, so the minimum is safe -- and never zero."""
    ky = abs(ZZ(ctx.discriminant))
    br = brunotte_norm_bound(ctx)
    bound = norm_bound_for(ctx)
    assert bound <= ky
    if br is not None and br > 0:
        assert bound == min(ky, br)
    assert bound >= 1


def test_every_result_is_within_the_norm_bound(ctx):
    bound = norm_bound_for(ctx)
    for y, _ in indecomposables_exhaustive(ctx):
        assert abs(ZZ(y.norm())) <= bound


@pytest.mark.slow
def test_ideal_sweep_is_ordered_by_norm(ctx):
    """
    Callers rely on increasing norm order: if x is decomposable then x = y + z
    with y indecomposable and N(y) < N(x), so everything needed to reject x has
    already been seen.
    """
    last = 0
    for nrm, _, _ in candidates_by_norm(ctx, min(norm_bound_for(ctx), 200)):
        assert nrm >= last
        last = nrm
