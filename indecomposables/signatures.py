"""
Canonical labelling of the signature classes of a totally real field.

Multiplication by a unit ``u`` is a bijection from the elements of signature
``sigma`` to those of signature ``sigma * sig(u)``, so sails and minimal
elements depend only on the coset of ``sigma`` modulo ``H = Sig(U)``.  The
classes are therefore ``G / H`` with ``G = F_2^n``.

Encoding
--------
A signature is a bitmask with ``+ -> 0``, ``- -> 1``, and ``sigma_1`` the **most
significant** bit.  The canonical representative of a class is the numerically
smallest mask in it, and classes are listed in ascending order of that mask.

Three consequences, all relied on downstream:

* class 0 is always the totally positive one, since ``0`` lies in every ``H``;
* every canonical representative begins with ``+``, because ``-1`` is a unit
  with signature ``(-,...,-)``, so ``sigma`` and ``-sigma`` always share a class
  and the one starting ``+`` has the smaller mask;
* the number of classes is ``2^(n - r)`` with ``r = unit_signature_rank``, which
  must also equal ``h^+ / h`` -- a cross-check over three independently computed
  columns.

Embedding order
---------------
Sign vectors only mean anything once ``sigma_1, ..., sigma_n`` are pinned down.
The convention, stated in ``data/README.md`` and assumed everywhere here, is
**increasing value of ``sigma_i(a)``** for ``a`` the root of the LMFDB defining
polynomial.  The labelling is canonical relative to the LMFDB label, in exactly
the way the ``basis`` column already is.

The group-theoretic half of this module is pure Python and testable without Sage;
only :func:`field_signature_classes` needs a field.
"""

from __future__ import annotations

__all__ = [
    "mask_to_vector", "vector_to_mask", "mask_to_string",
    "subgroup_closure", "signature_classes", "class_of",
    "field_signature_classes",
]


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------

def mask_to_vector(mask: int, n: int) -> tuple:
    """Bitmask to a vector of +-1, most significant bit first."""
    if not 0 <= mask < 2 ** n:
        raise ValueError(f"mask {mask} out of range for degree {n}")
    return tuple(-1 if (mask >> (n - 1 - i)) & 1 else 1 for i in range(n))


def vector_to_mask(vec) -> int:
    """Vector of signs (any nonzero reals, or +-1) to a bitmask."""
    mask = 0
    n = len(vec)
    for i, s in enumerate(vec):
        if s == 0:
            raise ValueError("zero has no signature")
        if s < 0:
            mask |= 1 << (n - 1 - i)
    return mask


def mask_to_string(mask: int, n: int) -> str:
    """Human-readable form, e.g. ``'+--'``.  For display only; not a column."""
    return "".join("+" if s > 0 else "-" for s in mask_to_vector(mask, n))


# ---------------------------------------------------------------------------
# Cosets
# ---------------------------------------------------------------------------

def subgroup_closure(generators, n: int) -> frozenset:
    """The F_2-subspace of ``F_2^n`` spanned by ``generators``, as masks."""
    H = {0}
    for g in generators:
        g = int(g)
        if g == 0 or g in H:
            continue
        H |= {h ^ g for h in H}
    return frozenset(H)


def signature_classes(unit_masks, n: int) -> tuple:
    """
    Canonical representatives of ``F_2^n / Sig(U)``, ascending.

    INPUT:

    - ``unit_masks`` -- signature masks of a generating set of the unit group,
      which must include ``-1`` (mask ``2^n - 1``)
    - ``n`` -- the degree

    OUTPUT: tuple of masks, ``classes[0] == 0``

    EXAMPLES::

        >>> signature_classes([0b11], 2)          # Sig(U) = {00, 11}
        (0, 1)
        >>> signature_classes([0b111], 3)         # only -1 is detected
        (0, 1, 2, 3)
        >>> signature_classes([0b111, 0b011], 3)
        (0, 1)
    """
    H = subgroup_closure(unit_masks, n)
    full = (1 << n) - 1
    if full not in H:
        raise ValueError(
            "the signature of -1 (all bits set) must lie in Sig(U); "
            "the generating set passed in is missing the torsion unit")
    seen, classes = set(), []
    for m in range(1 << n):                 # ascending, so the first hit is minimal
        if m in seen:
            continue
        seen |= {m ^ h for h in H}
        classes.append(m)
    expected = 1 << (n - _rank(H, n))
    if len(classes) != expected:
        raise AssertionError(f"got {len(classes)} classes, expected {expected}")
    return tuple(classes)


def _rank(H: frozenset, n: int) -> int:
    """Dimension of the subspace ``H`` over F_2, by Gaussian elimination."""
    basis, rank = [], 0
    for v in sorted(H, reverse=True):
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            rank += 1
    return rank


def class_of(mask: int, classes, n: int, unit_masks=None, H=None) -> int:
    """Index into ``classes`` of the class containing ``mask``."""
    if H is None:
        H = subgroup_closure(unit_masks or [], n)
    rep = min(mask ^ h for h in H)
    return classes.index(rep)


# ---------------------------------------------------------------------------
# Sage-facing
# ---------------------------------------------------------------------------

def field_signature_classes(ctx):
    """
    Canonical signature classes of the field behind ``ctx``.

    Returns ``(classes, vectors)`` where ``classes`` is the tuple of masks and
    ``vectors`` the matching +-1 vectors, suitable for the ``signature_masks``
    and ``signatures`` columns.
    """
    n = ctx.degree
    gens = [ctx.K(-1)] + list(ctx.fundamental_units)
    masks = [vector_to_mask(ctx.signature(u)) for u in gens]
    classes = signature_classes(masks, n)
    return classes, tuple(mask_to_vector(m, n) for m in classes)
