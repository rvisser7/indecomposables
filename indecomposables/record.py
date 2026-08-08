"""
The record layer, and the contract ``run_parallel.py`` depends on.

Three things live here: an **algorithm registry**, replacing the if-chain in the
old ``compute_indecomposables``; ``compute_record`` / ``failure_record``, which
produce immutable mappings satisfying :mod:`indecomposables.schema`; and thin
``to_row`` / ``parse_row`` wrappers so the driver never touches column names.

The record is a frozen mapping rather than a lazily-mutating object.  That is
what makes ``parse_row(to_row(r)) == r`` meaningful, and what reduces a worker's
resumable state to "which labels are done".

Partial results are first-class.  If some signature classes finish and others
blow up, the row is written with ``status = "partial"`` and the classes that did
finish are kept: throwing away good work because a later class failed would be
worse than recording exactly how far we got.
"""

from __future__ import annotations

from types import MappingProxyType

from . import schema
from .logging_utils import get_logger

logger = get_logger("record")

__all__ = [
    "register", "applicable_algorithms", "best_algorithm",
    "build_context", "compute_record", "failure_line", "to_row", "parse_row",
]


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY: list = []


def register(name, priority=0, complete=True):
    """
    Declare an algorithm.

    The decorated function has signature ``applies(ctx) -> callable | None``:
    return the routine to use, or ``None`` if it does not apply here.

    ``priority`` -- higher wins when several apply.  ``complete`` -- whether the
    routine returns *all* indecomposables.  A sail walk is sound but incomplete
    for degree at least 3, so it must never be selected as primary; keeping the
    flag on the registration rather than in a caller's head is what enforces it.
    """
    def deco(applies):
        _REGISTRY.append({"name": name, "priority": priority,
                          "complete": complete, "applies": applies,
                          "unavailable": None})
        _REGISTRY.sort(key=lambda e: -e["priority"])
        return applies
    return deco


def applicable_algorithms(ctx, complete_only=False):
    """
    Everything that applies, best first.  Cross-validation uses the whole list.

    An entry whose module will not import is skipped rather than allowed to
    propagate: a specialised family being unavailable must not stop a field from
    being computed by brute force.  Losing a fast path is a performance problem,
    not a correctness one -- the fallback is the exhaustive method.

    Unavailability is reported **once per algorithm per process**, not once per
    call.  It is a fixed fact about the installation, so repeating it for every
    field and every signature class would bury the actual results.
    """
    out = []
    for entry in _REGISTRY:
        if complete_only and not entry["complete"]:
            continue
        if entry["unavailable"] is not None:
            continue
        try:
            fn = entry["applies"](ctx)
        except (ImportError, NotImplementedError, AttributeError) as exc:
            entry["unavailable"] = str(exc) or type(exc).__name__
            logger.warning(
                "algorithm %r is not available and will be skipped for the rest "
                "of this run; falling back to the exhaustive method. (%s)",
                entry["name"], entry["unavailable"])
            continue
        if fn is not None:
            out.append((entry["name"], fn))
    return out


def unavailable_algorithms():
    """``{name: reason}`` for algorithms found unavailable during this run."""
    return {e["name"]: e["unavailable"] for e in _REGISTRY if e["unavailable"]}


def best_algorithm(ctx):
    found = applicable_algorithms(ctx, complete_only=True)
    if not found:
        raise LookupError(f"no complete algorithm applies to {ctx}")
    return found[0]




# ---------------------------------------------------------------------------
# Driver contract
# ---------------------------------------------------------------------------

def build_context(coeffs, label, known=None, **kw):
    from .context import FieldContext
    return FieldContext.from_coeffs(coeffs, label=label, known=known, **kw)


def _base(ctx):
    """Identity, provenance and joined invariants: everything not per-signature."""
    rec = {
        "lmfdb_label": ctx.label,
        "degree": int(ctx.degree),
        "discriminant": int(ctx.discriminant),
        "coeffs": [int(c) for c in ctx.K.defining_polynomial().list()[:-1]],
    }
    for name in ("regulator", "class_number", "narrow_class_number",
                 "is_monogenic", "monogenic_index", "num_subfields",
                 "is_galois", "galois_label", "family", "family_parameter"):
        v = ctx.known.get(name)
        if v is not None:
            rec[name] = v
    return rec


def compute_record(ctx, verify=True):
    """Compute one field across every signature class; return a frozen mapping."""
    from sage.all import ZZ

    best_algorithm(ctx)          # raises early if nothing applies
    rec = _base(ctx)

    masks, sigs, results, failures = _all_classes(ctx, verify)

    coords = [[[int(t) for t in ctx.coordinates(y)] for y, _, _ in cls]
              for cls in results]
    norms = [[int(abs(ZZ(ctx.K(y).norm()))) for y, _, _ in cls] for cls in results]
    sail = [[1 if on else 0 for _, on, _ in cls] for cls in results]

    traces = [[int(ZZ(ctx.K(y).trace())) for y, _, _ in cls] for cls in results]
    nonunit = [[n for n in cls if n != 1] for cls in norms]

    rec.update({
        "unit_signature_rank": int(ctx.unit_signature_rank),
        "num_signature_classes": len(sigs),
        "totally_positive_unit_index": len(sigs),
        "fundamental_units": [[int(t) for t in ctx.coordinates(u)]
                              for u in ctx.fundamental_units],
        "positive_unit_basis": [[int(t) for t in ctx.coordinates(u)]
                                for u in ctx.totally_positive_unit_basis],
        "signature_classes": [[int(s) for s in sig] for sig in sigs],

        # Everything below is one entry per signature class, in the order given
        # by signature_classes, with index 0 the totally positive one.
        "num_indecomposables": [len(c) for c in coords],
        "indecomposables": coords,
        "indecomposables_norms": norms,
        "indecomposables_traces": traces,
        "indecomposables_on_sail": sail,
        "min_norm_indecomposable": [min(c) if c else None for c in nonunit],
        "max_norm_indecomposable": [max(c) if c else None for c in norms],
        "min_trace_indecomposable": [min(c) if c else None for c in traces],
        "max_trace_indecomposable": [max(c) if c else None for c in traces],
        "num_on_sail": [sum(c) for c in sail],
    })
    _add_sail_columns(ctx, rec, sigs, results)
    return MappingProxyType(rec)


def _all_classes(ctx, verify):
    """
    Run every signature class, keeping whatever succeeds.

    A class that blows up is recorded in ``failures`` and left empty rather than
    aborting the field: a row covering three of four classes is worth more than
    no row at all, and ``status = "partial"`` says exactly that.
    """
    from .enumerate import indecomposables_in_class
    from .signatures import field_signature_classes

    masks, sigs = field_signature_classes(ctx)
    results, failures = [], []
    for sig in sigs:
        try:
            results.append(indecomposables_in_class(ctx, sig, verify=verify))
        except Exception as exc:                       # noqa: BLE001 - recorded
            logger.warning("%s signature %s failed: %s", ctx.label, sig, exc)
            failures.append(f"signature {tuple(sig)}: {type(exc).__name__}")
            results.append([])
    return masks, sigs, results, failures


def _add_sail_columns(ctx, rec, sigs, results):
    """
    Facet data, when :mod:`indecomposables.sail` is available.

    Vertices are stored as **indices into** ``indecomposables_all`` rather than
    repeated coordinates: half the size, and the sail columns then cannot
    disagree with the element columns about what a vertex is.
    """
    try:
        from .sail import sail_facets
    except (ImportError, NotImplementedError):
        return
    normals, vertices, counts, complete = [], [], [], []
    for sig, cls in zip(sigs, results):
        index = {tuple(ctx.coordinates(y)): i for i, (y, _, _) in enumerate(cls)}
        try:
            facets, ok = sail_facets(ctx, sig)
        except Exception as exc:                        # noqa: BLE001
            logger.warning("%s sail for %s failed: %s", ctx.label, sig, exc)
            normals.append([])
            vertices.append([])
            counts.append(0)
            complete.append(0)
            continue
        normals.append([[int(t) for t in ctx.coordinates(lam)] for lam, _ in facets])
        vertices.append([sorted(index[tuple(ctx.coordinates(v))] for v in verts)
                         for _, verts in facets])
        counts.append(len(facets))
        complete.append(1 if ok else 0)
    rec.update({
        "sail_facet_normals": normals,
        "sail_facet_vertices": vertices,
        "num_sail_facets_by_signature": counts,
        "sail_complete_by_signature": complete,
        "num_sail_facets": counts[0] if counts else None,
        "sail_complete": bool(complete[0]) if complete else None,
    })


def failure_line(label, reason):
    """
    One line of ``work/degree<n>/failures.txt`` for a field that did not finish.

    Failures are deliberately *not* rows.  Keeping them out means every column in
    the data file is an invariant of the field, and no column needs a NULL that
    means "not computed".  Nothing is lost: rows plus failures equals the scope,
    so a gap is still fully accounted for.
    """
    reason = " ".join(str(reason or "unspecified").split())[:200]
    return f"{label}|{reason}"


def to_row(record, tier=schema.FULL):
    return schema.encode_row(dict(record), tier)


def parse_row(line, tier=schema.FULL):
    return MappingProxyType(schema.decode_row(line, tier))


# ---------------------------------------------------------------------------
# Built-in registrations
# ---------------------------------------------------------------------------

@register("dress_scharlau", priority=100, complete=True)
def _dress_scharlau(ctx):
    if ctx.degree != 2:
        return None
    from .families.real_quadratic import indecomposables_dress_scharlau
    return indecomposables_dress_scharlau


@register("kala_tinkova", priority=90, complete=True)
def _kala_tinkova(ctx):
    """Simplest cubics; the parameter is read off the field, never used to build it."""
    if ctx.degree != 3:
        return None
    from .families.simplest_cubic import simplest_cubic_parameter
    if simplest_cubic_parameter(ctx) is None:
        return None
    from .families.simplest_cubic import indecomposables_kala_tinkova
    return indecomposables_kala_tinkova


@register("sail_walk", priority=50, complete=False)
def _sail_walk(ctx):
    """Sound but incomplete for degree >= 3, so never primary."""
    from .sail import sail_lattice_points
    return sail_lattice_points


@register("brute_force", priority=0, complete=True)
def _brute_force(ctx):
    from .enumerate import indecomposables_exhaustive
    return indecomposables_exhaustive
