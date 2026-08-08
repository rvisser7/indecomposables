"""
Finding all s-indecomposables, in every signature class.

Strategy
--------
The set of s-indecomposables is the **Hilbert basis** of the additive semigroup
``Lambda ∩ C``: the elements that are not sums of two others.  The cone ``C`` is
irrational with respect to ``Lambda``, so that basis is infinite -- but finite
modulo the totally positive units, which is what makes the computation finite.

Concretely, three steps.

1. **A rigorous cap.**  Kala-Yatsyna give ``N(alpha) <= |disc|`` for every
   indecomposable, and the minimal-trace representative of an orbit is balanced
   to within the covering radius ``rho`` of the log lattice of ``U^+``.  So each
   ``|sigma_i(alpha)| <= |disc|^(1/n) e^rho``, hence

       T2(alpha) = sum_i sigma_i(alpha)^2 <= n |disc|^(2/n) e^(2 rho).

   Unlike the trace, ``T2`` is positive definite on every signature class, so
   one bound covers them all.

2. **One lattice enumeration.**  ``T2`` is the trace form ``Tr(x y)``, integral
   and positive definite for a totally real field, so the candidates come from a
   single ``qfminim`` -- delegated to PARI rather than looped in Python, since
   this is the hot path.

3. **Tiered filtering**, cheapest first, in :mod:`indecomposables.certify`:
   the norm bound, then the codifferent certificate (which also settles sail
   membership), then the box test with an early exit.

A note on the module name: it shadows the Python builtin, so import the
submodule directly -- ``from indecomposables.enumerate import indecomposables_in_class``
-- and never ``from indecomposables import enumerate``, which would rebind
``enumerate`` for the rest of the importing file.  Every import inside the
package is relative and therefore unaffected.

The candidate count is roughly ``V_n n^(n/2) |disc|^(1/2) e^(n rho) / 2^n``.  The
``e^(n rho)`` is what limits the reachable degree: comfortable at 3, workable at
4 and 5, and only for small discriminant beyond that.  That is a property of the
problem, not of this implementation -- no efficient complete algorithm is known
for degree at least 3.
"""

from __future__ import annotations

from sage.all import Matrix, RealField, ZZ

from .certify import classify, is_indecomposable, is_on_sail
from .logging_utils import get_logger
from .signatures import field_signature_classes, vector_to_mask

logger = get_logger("enumerate")

__all__ = [
    "t2_bound", "trace_form", "candidates", "indecomposables_in_class",
    "indecomposables_exhaustive", "all_signature_classes",
]

#: Give up rather than enumerate more than this many candidates for one field.
MAX_CANDIDATES = 5_000_000


class CandidateExplosion(RuntimeError):
    """The candidate set exceeded :data:`MAX_CANDIDATES` for this field."""


# ---------------------------------------------------------------------------
# Bounds and the trace form
# ---------------------------------------------------------------------------

def t2_bound(ctx, norm_bound=None):
    """
    Cap on ``T2(x) = sum_i sigma_i(x)^2`` for a canonical s-indecomposable.

    Positive definite on every signature class, so a single bound serves all of
    them -- unlike the trace, which is only useful in the totally positive one.
    """
    R = RealField(ctx.prec)
    n = ctx.degree
    bound = R(abs(norm_bound if norm_bound is not None else ctx.discriminant))
    return R(n) * bound ** (R(2) / n) * (2 * ctx.log_unit_covering_radius).exp()


def trace_form(ctx):
    """Gram matrix ``Tr(b_i b_j)`` of the integral basis: integral and positive definite."""
    if getattr(ctx, "_trace_form", None) is None:
        b = ctx.basis
        ctx._trace_form = Matrix(ZZ, [[ZZ((ctx.K(bi) * ctx.K(bj)).trace())
                                       for bj in b] for bi in b])
    return ctx._trace_form


def candidates(ctx, bound=None, cap=MAX_CANDIDATES):
    """
    Every nonzero ``x`` in the order with ``T2(x) <= bound``, up to sign pairs.

    ``qfminim`` returns one of each ``+-`` pair, so callers that care about
    signature classes must consider both ``x`` and ``-x``.  Delegated to PARI:
    this is the inner loop and must not be a Python loop.
    """
    G = trace_form(ctx)
    B = bound if bound is not None else t2_bound(ctx)
    count, _, vectors = G.__pari__().qfminim(int(B) + 1, cap, flag=0)
    if int(count) >= cap:
        raise CandidateExplosion(
            f"{ctx.label or ctx.K}: more than {cap} candidates with T2 <= {float(B):.3g}; "
            "the regulator makes exhaustive search impractical for this field")
    logger.info("%s: %s candidate pairs with T2 <= %.4g",
                ctx.label or ctx.K, count, float(B))
    for col in Matrix(ZZ, vectors.mattranspose().sage()).rows():
        yield ctx.from_coordinates(col)


# ---------------------------------------------------------------------------
# Minimal elements
# ---------------------------------------------------------------------------

def indecomposables_in_class(ctx, signature=None, bound=None, verify=True):
    """
    The s-indecomposables of one signature class, canonical and deduplicated.

    Returns a list of ``(element, on_sail, reason)``.  ``reason`` records which
    test decided the element, since only ``box_exhausted`` depends on a
    numerical search having been complete.
    """
    n = ctx.degree
    if signature is None:
        signature = [1] * n
    signature = tuple(int(s) for s in signature)
    totally_positive = all(s > 0 for s in signature)
    norm_bound = abs(ZZ(bound if bound is not None else ctx.discriminant))
    embs = ctx.real_embeddings()

    found = []
    for x in candidates(ctx, t2_bound(ctx, norm_bound)):
        for y in (x, -x):
            vals = [embs[i](y) for i in range(n)]
            if any(v == 0 for v in vals):
                continue
            if tuple(1 if v > 0 else -1 for v in vals) != signature:
                continue
            if abs(ZZ(y.norm())) > norm_bound:
                continue
            if totally_positive:
                indec, on_sail, reason, _ = classify(y, ctx, verify=verify)
                if indec:
                    found.append((y, on_sail, reason))
            elif is_indecomposable(y, ctx, signature, verify=verify):
                found.append((y, is_on_sail_general(y, ctx), "box_exhausted"))

    canon = ctx.normalizer
    seen = {}
    for y, on_sail, reason in found:
        c = canon.canonical(y) if totally_positive else y
        key = tuple(ctx.coordinates(c))
        if key not in seen or on_sail:
            seen[key] = (c, on_sail, reason)
    out = sorted(seen.values(), key=lambda t: canon.sort_key(t[0]))
    logger.info("%s signature %s: %s s-indecomposable, %s on the sail",
                ctx.label or ctx.K, signature, len(out), sum(1 for t in out if t[1]))
    return out


def is_on_sail_general(x, ctx):
    """
    Sail membership in a general signature class.

    The certificate argument is stated for the totally positive cone; for other
    classes the supporting hyperplane is found by :mod:`indecomposables.sail`
    from the polytope, so this is a placeholder returning ``False`` rather than
    a wrong ``True``.  Being conservative here means ``num_on_sail`` can only
    ever undercount, never claim a point is on the sail when it is not.
    """
    if all(ctx.real_embeddings()[i](x) > 0 for i in range(ctx.degree)):
        return is_on_sail(x, ctx)
    return False


# ---------------------------------------------------------------------------
# Entry points used by the registry
# ---------------------------------------------------------------------------

def indecomposables_exhaustive(ctx, verify=True):
    """All indecomposables, i.e. the s-indecomposables of the totally positive class."""
    return [(y, on_sail) for y, on_sail, _ in indecomposables_in_class(ctx, verify=verify)]


def all_signature_classes(ctx, verify=True):
    """
    The s-indecomposables of every signature class.

    Returns ``(masks, signatures, results)`` with ``results`` a list parallel to
    ``signatures``, index 0 being the totally positive class -- matching the
    ordering the ``_by_signature`` schema columns assume.
    """
    masks, sigs = field_signature_classes(ctx)
    results = []
    for sig in sigs:
        results.append(indecomposables_in_class(ctx, sig, verify=verify))
    assert masks[0] == 0 and vector_to_mask(sigs[0]) == 0, \
        "class 0 must be the totally positive one"
    return masks, sigs, results
