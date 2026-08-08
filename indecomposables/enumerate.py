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
    "t2_bound", "trace_form", "candidates", "candidates_by_norm",
    "brunotte_norm_bound", "norm_bound_for", "indecomposables_in_class",
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
    Coefficient vectors of every nonzero ``x`` with ``T2(x) <= bound``.

    Yields **integer coordinate vectors**, not field elements.  Constructing a
    Sage element costs far more than the arithmetic needed to reject it, and the
    overwhelming majority are rejected: for Q(sqrt 61) the bound admits some
    80,000 vectors of which about 74 survive the norm filter.  So the caller
    screens in coordinate space and builds elements only for survivors.

    ``qfminim`` returns one of each ``+-`` pair, so a caller that cares about
    signature classes must consider both ``v`` and ``-v``.
    """
    G = trace_form(ctx)
    B = bound if bound is not None else t2_bound(ctx)
    count, _, vectors = G.__pari__().qfminim(int(B) + 1, cap, flag=0)
    if int(count) >= cap:
        raise CandidateExplosion(
            f"{ctx.label or ctx.K}: more than {cap} candidates with "
            f"T2 <= {float(B):.3g}; the regulator makes exhaustive search "
            "impractical for this field")
    logger.info("%s: %s candidate pairs with T2 <= %.4g",
                ctx.label or ctx.K, count, float(B))
    return Matrix(ZZ, vectors.mattranspose().sage()).rows()


def unit_representatives(ctx):
    """
    Coset representatives of ``U / U^2``: the ``2^n`` products of ``-1`` and the
    fundamental units with exponents in ``{0, 1}``.

    Multiplying a generator by these reaches every signature its associates can
    have, which is what lets one sweep over ideals cover all signature classes.
    """
    if getattr(ctx, "_unit_reps", None) is None:
        from sage.all import prod
        gens = [ctx.K(-1)] + list(ctx.fundamental_units)
        reps = []
        for mask in range(2 ** len(gens)):
            reps.append(prod((g for i, g in enumerate(gens) if mask >> i & 1),
                             ctx.K.one()))
        ctx._unit_reps = reps
    return ctx._unit_reps


def candidates_by_norm(ctx, norm_bound, chunk=256):
    """
    Totally positive generators of principal ideals of norm at most ``norm_bound``.

    Yields ``(norm, element, signature)`` for every associate of every principal
    ideal generator, one per signature the associate can have.

    This is the good way to generate candidates, and it is why the regulator
    does not appear.  A bounded region of Minkowski space has to contain a
    fundamental domain for the unit action, so :func:`candidates` pays a factor
    ``e^(n rho)`` -- for Q(sqrt 61) that is 74,738 lattice vectors against 61
    ideal norms.  An **ideal is already unit-invariant**, so enumerating by norm
    gives one object per orbit for free.

    Norms are processed in increasing order, which callers rely on: if ``x`` is
    decomposable then ``x = y + z`` with ``y`` indecomposable and ``N(y) <
    N(x)``, so everything needed to reject ``x`` has already been seen.
    """
    from sage.rings.number_field.bdd_height import bdd_norm_pr_ideal_gens

    embs = ctx.real_embeddings()
    reps = unit_representatives(ctx)
    n = ctx.degree

    for lo in range(1, int(norm_bound) + 1, chunk):
        hi = min(lo + chunk - 1, int(norm_bound))
        gens = bdd_norm_pr_ideal_gens(ctx.K, range(lo, hi + 1))
        for nrm in range(lo, hi + 1):
            for g in gens.get(nrm, []):
                seen = set()
                for u in reps:
                    y = ctx.K(g) * u
                    vals = [embs[i](y) for i in range(n)]
                    if any(v == 0 for v in vals):
                        continue
                    sig = tuple(1 if v > 0 else -1 for v in vals)
                    if sig in seen:
                        continue
                    seen.add(sig)
                    yield ZZ(nrm), y, sig


def brunotte_norm_bound(ctx):
    """
    Brunotte's bound on the norm of an indecomposable.

    ``max_i prod_rho max(1, |sigma_i(u_rho)|)`` raised to the degree, over a
    basis of the totally positive units.  Sometimes sharper than Kala-Yatsyna's
    ``|disc|`` and sometimes far worse, so callers take the minimum of the two.

    Brunotte, *Zur Zerlegung totalpositiver Zahlen in Ordnungen totalreeller
    algebraischer Zahlkoerper*, Arch. Math. 41 (1983).
    """
    from sage.all import RealField
    R = RealField(ctx.prec)
    embs = ctx.real_embeddings()
    units = ctx.totally_positive_unit_basis
    if not units:
        return None
    try:
        c = max(prod_over(max(R(1), abs(R(embs[i](u)))) for u in units)
                for i in range(ctx.degree))
        return ZZ((c ** ctx.degree).ceil())
    except (ValueError, ZeroDivisionError):
        return None


def prod_over(values):
    out = None
    for v in values:
        out = v if out is None else out * v
    return out


def _screen(ctx, norm_bound, signature):
    """
    A fast approximate filter on coordinate vectors.

    Returns a callable ``(vector, sign) -> bool``: True means "might qualify,
    build the element and check exactly".  Everything here is double precision
    on a precomputed embedding matrix -- no Sage elements, no high-precision
    arithmetic.

    Soundness comes from tolerances that always err towards *keeping* a vector:
    a value too close to zero to sign reliably, or a norm too close to the
    bound, is passed through to the exact test.  So the screen can waste work
    but cannot discard a genuine s-indecomposable.
    """
    n = ctx.degree
    embs = ctx.real_embeddings(53)
    M = [[float(e(b)) for e in embs] for b in ctx.basis]
    scale = max(abs(x) for row in M for x in row) or 1.0
    zero_tol = 1e-9 * scale
    slack = 1.0 + 1e-9
    want = tuple(int(s) for s in signature)

    def keep(v, sign):
        vals = [sign * sum(float(c) * M[j][i] for j, c in enumerate(v))
                for i in range(n)]
        prod = 1.0
        for i, x in enumerate(vals):
            if abs(x) <= zero_tol:
                return True                 # cannot sign it; let exact decide
            if (x > 0) != (want[i] > 0):
                return False
            prod *= abs(x)
        return prod <= norm_bound * slack

    return keep


# ---------------------------------------------------------------------------
# Minimal elements
# ---------------------------------------------------------------------------

def norm_bound_for(ctx, bound=None):
    """
    Norm bound for an s-indecomposable: the best of the two known bounds.

    Kala-Yatsyna give ``|disc|``; Brunotte's is sometimes sharper.  Both are
    proved, so taking the minimum is safe.
    """
    if bound is not None:
        return abs(ZZ(bound))
    ky = abs(ZZ(ctx.discriminant))
    br = brunotte_norm_bound(ctx)
    return min(ky, br) if br is not None and br > 0 else ky


def indecomposables_in_class(ctx, signature=None, bound=None, verify=True,
                             sail=True, method="ideals"):
    """
    The s-indecomposables of one signature class, canonical and deduplicated.

    Returns a list of ``(element, on_sail, reason)``.  ``reason`` records which
    test decided the element, since only ``box_exhausted`` depends on a
    numerical search having been complete.

    With ``sail=False`` no sail membership is determined and ``on_sail`` is
    ``None`` throughout, so every sail column ends up null.

    ``method`` selects candidate generation:

    ``"ideals"``
        Principal ideals of bounded norm.  The regulator does not appear, so
        this is the one to use -- see :func:`candidates_by_norm`.
    ``"lattice"``
        Lattice points of bounded ``T2``.  Costs a factor ``e^(n rho)`` more,
        but is completely independent of the ideal machinery, which makes it a
        useful second opinion: the two must find the same set.
    """
    n = ctx.degree
    if signature is None:
        signature = [1] * n
    signature = tuple(int(s) for s in signature)
    totally_positive = all(s > 0 for s in signature)
    norm_bound = norm_bound_for(ctx, bound)

    if method == "ideals":
        found = _by_ideals(ctx, signature, norm_bound, totally_positive, verify, sail)
    elif method == "lattice":
        found = _by_lattice(ctx, signature, norm_bound, totally_positive, verify, sail)
    else:
        raise ValueError(f"unknown method {method!r}; use 'ideals' or 'lattice'")

    canon = ctx.normalizer
    seen = {}
    for y, on_sail, reason in found:
        c = canon.canonical(y) if totally_positive else y
        key = tuple(ctx.coordinates(c))
        if key not in seen or on_sail:
            seen[key] = (c, on_sail, reason)
    out = sorted(seen.values(), key=lambda t: canon.sort_key(t[0]))
    logger.info("%s signature %s: %s s-indecomposable%s (%s, norm bound %s)",
                ctx.label or ctx.K, signature, len(out),
                f", {sum(1 for t in out if t[1])} on the sail" if sail else
                " (sail not computed)", method, norm_bound)
    return out


def _decide(y, ctx, signature, totally_positive, verify, sail):
    """Classify one candidate; returns ``(element, on_sail, reason)`` or ``None``."""
    if totally_positive and sail:
        indec, on_sail, reason, _ = classify(y, ctx, verify=verify)
        return (y, on_sail, reason) if indec else None
    if is_indecomposable(y, ctx, signature, verify=verify):
        return (y, is_on_sail_general(y, ctx) if sail else None, "box_exhausted")
    return None


def _by_ideals(ctx, signature, norm_bound, totally_positive, verify, sail):
    """
    Candidates from principal ideals of bounded norm.

    Also records, per signature, the smallest absolute norm attained -- which is
    what ``min_norm_all_signatures`` and ``min_norm_fprimea_signature`` are, and
    it falls out of this sweep for nothing.  The lattice method cannot supply it
    as cheaply, so it is only populated here.
    """
    found = []
    smallest = getattr(ctx, "_smallest_norm_by_signature", None)
    if smallest is None:
        smallest = ctx._smallest_norm_by_signature = {}
    for nrm, y, sig in candidates_by_norm(ctx, norm_bound):
        prev = smallest.get(sig)
        if prev is None or nrm < prev:
            smallest[sig] = nrm
        if sig != signature:
            continue
        decided = _decide(y, ctx, signature, totally_positive, verify, sail)
        if decided is not None:
            found.append(decided)
    return found


def _by_lattice(ctx, signature, norm_bound, totally_positive, verify, sail):
    """Candidates from lattice points of bounded T2 -- the independent check."""
    n = ctx.degree
    embs = ctx.real_embeddings()
    keep = _screen(ctx, norm_bound, signature)
    found, screened, exact = [], 0, 0
    for v in candidates(ctx, t2_bound(ctx, norm_bound)):
        for sign in (1, -1):
            screened += 1
            if not keep(v, sign):
                continue
            exact += 1
            y = ctx.from_coordinates([sign * c for c in v])
            vals = [embs[i](y) for i in range(n)]
            if any(v_ == 0 for v_ in vals):
                continue
            if tuple(1 if v_ > 0 else -1 for v_ in vals) != signature:
                continue
            if abs(ZZ(y.norm())) > norm_bound:
                continue
            decided = _decide(y, ctx, signature, totally_positive, verify, sail)
            if decided is not None:
                found.append(decided)
    logger.info("%s signature %s: screened %s, examined %s exactly (%.2f%%)",
                ctx.label or ctx.K, signature, screened, exact,
                100.0 * exact / max(screened, 1))
    return found


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

def indecomposables_exhaustive(ctx, verify=True, sail=True, method="ideals"):
    """All indecomposables, i.e. the s-indecomposables of the totally positive class."""
    return [(y, on_sail)
            for y, on_sail, _ in indecomposables_in_class(
                ctx, verify=verify, sail=sail, method=method)]


def all_signature_classes(ctx, verify=True, sail=True, method="ideals"):
    """
    The s-indecomposables of every signature class.

    Returns ``(masks, signatures, results)`` with ``results`` a list parallel to
    ``signatures``, index 0 being the totally positive class -- matching the
    ordering the ``_by_signature`` schema columns assume.
    """
    masks, sigs = field_signature_classes(ctx)
    results = []
    for sig in sigs:
        results.append(indecomposables_in_class(ctx, sig, verify=verify,
                                                sail=sail, method=method))
    assert masks[0] == 0 and vector_to_mask(sigs[0]) == 0, \
        "class 0 must be the totally positive one"
    return masks, sigs, results
