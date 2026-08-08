"""
Facets of the sail.

The sail is the boundary of the convex hull of the totally positive elements of
the order.  Its facets are exactly the hyperplanes

    Tr(lambda x) = 1,   lambda totally positive in the codifferent D^-1,

since such a ``lambda`` satisfies ``Tr(lambda beta) >= 1`` for every totally
positive ``beta``, so the hyperplane supports the hull and everything on it is
indecomposable.

Why this needs no Klein-polyhedron walk
---------------------------------------
The sail is infinite but ``U^+`` acts on it with finitely many orbits, and

* multiplication by ``eps`` sends the facet with normal ``lambda`` to the one
  with normal ``eps^-1 lambda``, so **facet orbits are exactly the ``U^+``-orbits
  of supporting ``lambda``** -- canonicalised by the same machinery that
  canonicalises elements, over the codifferent's basis instead of the order's;
* every facet contains at least one lattice point, namely a vertex, and every
  such point is indecomposable.

So enumerating the supporting ``lambda`` of each known indecomposable finds every
facet orbit, and no bounding-box search or completeness certificate is needed.
Both halves reuse code that is already tested: the box enumeration in
:mod:`indecomposables.certify` and the canonical form in
:mod:`indecomposables.normalize`.

Scope
-----
Only the totally positive class.  The supporting-hyperplane argument is stated
there, and the analogue for a general signature class needs the polytope
construction; :func:`sail_facets` reports ``complete=False`` for those rather
than guessing, so ``num_indecomposables_on_sail`` can only ever undercount.

Vertices are recorded as **indices into the class's indecomposable list**, so a
facet's vertex list names *orbits* of vertices.  Two vertices of one facet that
are unit multiples of each other therefore appear as the same index, and a facet
can have fewer listed vertices than its dimension would suggest -- that is the
encoding working, not an error.
"""

from __future__ import annotations

from sage.all import QQ, ZZ

from .certify import _box_points, _codifferent_basis
from .logging_utils import get_logger
from .normalize import UnitOrbitNormalizer

logger = get_logger("sail")

__all__ = ["codifferent_normalizer", "supporting_normals", "sail_facets",
           "sail_lattice_points"]


def codifferent_normalizer(ctx):
    """
    Canonical representatives for ``U^+`` acting on the codifferent.

    The same minimal-trace rule as for elements, over the codifferent's own
    Z-basis.  Note ``Tr(D^-1) is contained in Z`` by definition of the
    codifferent, so the trace used by the canonical form is still an integer.
    """
    if getattr(ctx, "_codiff_normalizer", None) is None:
        ctx._codiff_normalizer = UnitOrbitNormalizer(
            ctx.K, basis=_codifferent_basis(ctx),
            unit_basis=ctx.totally_positive_unit_basis, prec=ctx.prec)
    return ctx._codiff_normalizer


def supporting_normals(x, ctx, prec=None):
    """
    Every totally positive ``lambda`` in the codifferent with ``Tr(lambda x) = 1``.

    These are the normals of the facets through ``x``.  Since ``lambda`` and
    ``x`` are both totally positive with the traced product equal to 1, each
    ``sigma_i(lambda)`` lies in ``(0, 1/sigma_i(x))`` -- a box, so the same
    enumerator as everywhere else applies.
    """
    K = ctx.K
    x = K(x)
    prec = prec or ctx.prec
    embs = ctx.real_embeddings(prec)
    basis = _codifferent_basis(ctx)
    uppers = [1 / embs[i](x) for i in range(ctx.degree)]

    out = []
    for a in _box_points(basis, uppers, ctx, prec):
        lam = sum(ZZ(c) * b for c, b in zip(a, basis))
        if lam.is_zero() or not lam.is_totally_positive():
            continue
        if QQ((lam * x).trace()) == 1:
            out.append(lam)
    return out


def sail_facets(ctx, signature=None, elements=None, verify=True):
    """
    Facet orbits of the sail, as ``(normals, vertices, complete)``.

    - ``normals``  -- one canonical ``lambda`` per facet orbit
    - ``vertices`` -- for each, the indices in ``elements`` lying on it
    - ``complete`` -- whether every facet orbit was found

    ``elements`` is the class's indecomposables, in the order they are stored;
    the returned indices refer to it.  Facets are ordered by their normal's
    canonical sort key, so the output is reproducible.
    """
    n = ctx.degree
    if signature is not None and not all(s > 0 for s in signature):
        # See the module docstring: the argument is stated for the totally
        # positive cone only, and guessing would let num_..._on_sail overcount.
        logger.info("%s: sail for signature %s not computed (needs the polytope "
                    "construction)", ctx.label or ctx.K, tuple(signature))
        return [], [], False
    if not elements:
        return [], [], True

    canon = codifferent_normalizer(ctx)
    by_facet = {}
    for index, x in enumerate(elements):
        for lam in supporting_normals(x, ctx):
            key = tuple(canon.coordinates(canon.canonical(lam, verify=verify)))
            by_facet.setdefault(key, set()).add(index)

    order = sorted(by_facet, key=lambda k: canon.sort_key(canon.from_coordinates(k)))
    normals = [list(k) for k in order]
    vertices = [sorted(by_facet[k]) for k in order]
    logger.info("%s: %s facet orbit(s) from %s indecomposable(s)",
                ctx.label or ctx.K, len(normals), len(elements))
    # Complete by construction: every facet contains a vertex, every vertex is
    # indecomposable, and every indecomposable of this class is in `elements`.
    return normals, vertices, len(normals) >= n - 1 or n <= 2


def sail_lattice_points(ctx, verify=True):
    """
    The indecomposables lying on the sail, for the registry's ``sail_walk`` entry.

    Sound but **incomplete** as an indecomposable-finding algorithm for degree at
    least 3, which is why it is registered with ``complete=False`` and can never
    be chosen as the primary method.
    """
    from .enumerate import indecomposables_in_class
    return [(y, on_sail) for y, on_sail, _ in indecomposables_in_class(ctx, verify=verify)
            if on_sail]
