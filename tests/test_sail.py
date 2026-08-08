"""
Sail facets.

The implementation rests on two claims, and these tests check both:

1. the facets of the sail are the hyperplanes ``Tr(lambda x) = 1`` for totally
   positive ``lambda`` in the codifferent, and facet orbits under ``U^+``
   correspond to ``U^+``-orbits of such ``lambda``;
2. enumerating those ``lambda`` over the known indecomposables finds every facet,
   because every facet contains a vertex and every vertex is indecomposable.

The degree-2 test below is the sharp one.  There the sail is a polygon, so the
number of edge orbits must equal the number of vertex orbits.  A mismatch would
mean claim 1 is wrong -- specifically that canonical ``lambda`` does not identify
facet orbits -- so this is the test that decides whether the approach is sound.
I could only check it in a float model, where orbit hashing and truncation made
the answer ambiguous; under exact arithmetic it is decisive.
"""

import pytest

sage = pytest.importorskip("sage.all", reason="Sage not available")
NumberField = sage.NumberField
PolynomialRing = sage.PolynomialRing
QQ, ZZ = sage.QQ, sage.ZZ

pytestmark = [pytest.mark.sage, pytest.mark.sail]

from indecomposables.context import FieldContext                    # noqa: E402
from indecomposables.certify import is_indecomposable, is_on_sail    # noqa: E402
from indecomposables.enumerate import indecomposables_in_class       # noqa: E402
from indecomposables.sail import (                                   # noqa: E402
    codifferent_normalizer, sail_facets, supporting_normals,
)

R = PolynomialRing(QQ, "x")
x = R.gen()

QUADRATIC = [5, 2, 3, 13, 6, 7, 10, 21, 33]
CUBIC = {
    "3.3.49.1": x ** 3 - x ** 2 - 2 * x + 1,
    "3.3.81.1": x ** 3 - 3 * x - 1,
    "3.3.148.1": x ** 3 - x ** 2 - 3 * x + 1,
}


def _quadratic(D):
    f = x ** 2 - x - (D - 1) // 4 if D % 4 == 1 else x ** 2 - D
    return FieldContext(NumberField(f, "a"), label=f"Q(sqrt {D})")


@pytest.fixture(scope="session", params=sorted(CUBIC))
def cubic(request):
    return FieldContext(NumberField(CUBIC[request.param], "a"), label=request.param)


# ---------------------------------------------------------------------------
# The supporting hyperplanes themselves
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", QUADRATIC)
def test_normals_are_totally_positive_codifferent_elements(D):
    ctx = _quadratic(D)
    for y, _, _ in indecomposables_in_class(ctx):
        for lam in supporting_normals(y, ctx):
            assert lam.is_totally_positive()
            assert QQ((lam * y).trace()) == 1
            # in the codifferent: Tr(lambda * O_K) is contained in Z
            assert all(QQ((lam * b).trace()) in ZZ for b in ctx.basis)


@pytest.mark.parametrize("D", QUADRATIC)
def test_a_normal_exists_exactly_when_on_the_sail(D):
    ctx = _quadratic(D)
    for y, _, _ in indecomposables_in_class(ctx):
        assert bool(supporting_normals(y, ctx)) == is_on_sail(y, ctx)


def test_every_vertex_of_every_facet_is_indecomposable(cubic):
    elements = [y for y, _, _ in indecomposables_in_class(cubic)]
    normals, vertices, _ = sail_facets(cubic, None, elements)
    for verts in vertices:
        for i in verts:
            assert is_indecomposable(elements[i], cubic)


# ---------------------------------------------------------------------------
# The claim the approach rests on
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", QUADRATIC)
def test_degree_two_facet_orbits_equal_vertex_orbits(D):
    """
    In degree 2 the sail is a polygon, so edges and vertices are in bijection
    and their orbit counts must agree.

    If this fails, canonical lambda is not identifying facet orbits and
    sail_facets is wrong -- not merely imprecise.
    """
    ctx = _quadratic(D)
    elements = [y for y, _, _ in indecomposables_in_class(ctx)]
    normals, vertices, _ = sail_facets(ctx, None, elements)
    on_sail = [y for y in elements if is_on_sail(y, ctx)]
    assert len(normals) == len(on_sail), (
        f"Q(sqrt {D}): {len(normals)} facet orbits but {len(on_sail)} sail "
        "vertex orbits; in degree 2 these must be equal")


@pytest.mark.parametrize("D", QUADRATIC)
def test_each_degree_two_facet_has_two_vertices_counted_with_orbits(D):
    """
    An edge has two endpoints, but they may lie in the same orbit, so the stored
    vertex list has one or two entries -- never more.
    """
    ctx = _quadratic(D)
    elements = [y for y, _, _ in indecomposables_in_class(ctx)]
    _, vertices, _ = sail_facets(ctx, None, elements)
    for verts in vertices:
        assert 1 <= len(verts) <= 2, verts


# ---------------------------------------------------------------------------
# Bookkeeping
# ---------------------------------------------------------------------------

def test_facets_are_deduplicated_and_ordered(cubic):
    elements = [y for y, _, _ in indecomposables_in_class(cubic)]
    normals, vertices, _ = sail_facets(cubic, None, elements)
    assert len(normals) == len({tuple(v) for v in normals})
    canon = codifferent_normalizer(cubic)
    keys = [canon.sort_key(canon.from_coordinates(v)) for v in normals]
    assert keys == sorted(keys)
    assert len(vertices) == len(normals)


def test_vertex_indices_are_in_range(cubic):
    elements = [y for y, _, _ in indecomposables_in_class(cubic)]
    _, vertices, _ = sail_facets(cubic, None, elements)
    for verts in vertices:
        assert verts == sorted(set(verts))
        assert all(0 <= i < len(elements) for i in verts)


def test_non_totally_positive_class_is_reported_incomplete(cubic):
    """
    The supporting-hyperplane argument is only stated for the totally positive
    cone, so other classes must report incomplete rather than guess.
    """
    from indecomposables.signatures import field_signature_classes
    _, sigs = field_signature_classes(cubic)
    for sig in sigs:
        if all(s > 0 for s in sig):
            continue
        normals, vertices, complete = sail_facets(cubic, sig, [])
        assert normals == [] and vertices == [] and complete is False


def test_on_sail_flags_agree_with_the_facets(cubic):
    """An element is flagged on-sail exactly when some facet contains it."""
    results = indecomposables_in_class(cubic)
    elements = [y for y, _, _ in results]
    flags = [on for _, on, _ in results]
    _, vertices, _ = sail_facets(cubic, None, elements)
    covered = {i for verts in vertices for i in verts}
    for i, flag in enumerate(flags):
        assert flag == (i in covered), (
            f"element {i} flagged on_sail={flag} but "
            f"{'is' if i in covered else 'is not'} on a facet")
