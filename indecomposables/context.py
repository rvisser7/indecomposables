"""
Per-field context: everything expensive that is shared across one field's work.

This is the half of the old ``NumberFieldData`` that genuinely deserves to be an
object.  It is built once per field and passed to the algorithm functions, which
are otherwise free functions taking Sage objects.

What it deliberately does *not* hold: the indecomposables, the output row, the
logger configuration, or a ``verbose`` flag.  Those belong to ``record.py`` and
to the caller respectively.  Fusing them is what made ``from_data_row``
unwritable in the original class.

Two rules that matter for correctness rather than tidiness:

1. **Signs are exact, logs are not.**  Total positivity and signature vectors go
   through ``K.embeddings(AA)``; logarithms and geometry go through
   ``K.real_embeddings(prec)``.  The original code used
   ``element.complex_embeddings(prec)`` for both -- ordering complex numbers
   against zero, and taking ``.log()`` on a complex value, where a spurious
   imaginary part silently selects the principal branch instead of raising.

2. **One field, one basis.**  The canonical form of an indecomposable depends on
   the integral basis, so a field must never be built two different ways within
   a run.  ``FieldContext.from_coeffs`` is the single door.
"""

from __future__ import annotations

from functools import cached_property

from sage.all import AA, Matrix, NumberField, PolynomialRing, QQ, RealField, ZZ

from .logging_utils import get_logger

ZZx = PolynomialRing(ZZ, "x")

logger = get_logger("context")


def _try(compute):
    """
    Run ``compute``, returning ``None`` if it fails.

    Used for the LMFDB-style invariants.  These are conveniences, not the point
    of the computation, so one of them being slow or unavailable must not lose
    the indecomposables for that field -- a NULL in the column is the right
    outcome, and it is recoverable by supplying the value in the input table.
    """
    try:
        return compute()
    except Exception as exc:                       # noqa: BLE001 - recorded
        logger.info("invariant unavailable: %s: %s", type(exc).__name__, exc)
        return None

__all__ = ["FieldContext"]


class FieldContext:
    """
    Cached invariants and geometry for one totally real field.

    INPUT:

    - ``K`` -- a totally real number field
    - ``label`` -- LMFDB label, or ``None``
    - ``order`` -- an order in ``K`` (default: the maximal order)
    - ``basis`` -- integral basis fixing coefficient vectors and the canonical
      form.  Pass LMFDB's ``zk`` in production so that stored rows are
      reproducible; the default is the order's own basis.
    - ``known`` -- invariants already known from LMFDB (``regulator``,
      ``class_number``, ``is_monogenic``, ``num_subfields``, ...).  Anything
      supplied here is *not* recomputed: recomputing regulators and class
      numbers across 200k fields is a large fraction of total runtime for values
      LMFDB already publishes.
    - ``prec`` -- working precision in bits for the inexact geometry
    """

    def __init__(self, K, label=None, order=None, basis=None, known=None,
                 prec=128, basis_name="order"):
        self.K = K
        self.label = label
        self.prec = prec
        self.known = dict(known or {})
        self.basis_name = basis_name

        self.embeddings_exact = K.embeddings(AA)
        if len(self.embeddings_exact) != K.degree():
            raise ValueError(f"{K} is not totally real")

        self.order = order if order is not None else K.maximal_order()
        self.basis = list(basis) if basis is not None else list(self.order.basis())
        if len(self.basis) != K.degree():
            raise ValueError("integral basis has the wrong length")

    # -- construction -------------------------------------------------------

    @classmethod
    def from_coeffs(cls, coeffs, label=None, **kw):
        """
        Build from a coefficient list, constant term first, leading 1 optional.

        The single construction path.  Family-specific code (simplest cubics,
        biquadratics) must recover its parameter *from* the field rather than
        building its own copy, or the two copies can disagree on the defining
        polynomial and hence on the canonical representatives.
        """
        coeffs = [ZZ(c) for c in coeffs]
        if coeffs[-1] != 1:
            coeffs = coeffs + [ZZ(1)]
        return cls(NumberField(ZZx(coeffs), names="a"), label=label, **kw)

    # -- cheap invariants ---------------------------------------------------

    @property
    def degree(self):
        return self.K.degree()

    @cached_property
    def discriminant(self):
        return ZZ(self.known.get("discriminant") or self.K.discriminant())

    @cached_property
    def regulator(self):
        v = self.known.get("regulator")
        return float(v) if v is not None else _try(
            lambda: float(self.K.regulator()))

    @cached_property
    def class_number(self):
        """
        Class number.

        Computed only when not supplied.  At higher degree this is one of the
        slowest things here, so joining it from LMFDB is worth doing for a large
        run -- see the note in totally_real_fields/README.md.
        """
        v = self.known.get("class_number")
        return ZZ(v) if v is not None else _try(lambda: ZZ(self.K.class_number()))

    @cached_property
    def narrow_class_number(self):
        """
        Narrow class number h+.

        Computed from the narrow class group rather than derived as
        ``h * 2^(n-r)``: deriving it would make the ``h+/h ==
        num_signature_classes`` consistency check tautological, and that check
        has already earned its place.
        """
        v = self.known.get("narrow_class_number")
        if v is not None:
            return ZZ(v)
        return _try(lambda: ZZ(self.K.narrow_class_group().order()))

    @cached_property
    def monogenic_index(self):
        """``[O_K : Z[a]]``, the index of the power order in the maximal order."""
        v = self.known.get("monogenic_index")
        if v is not None:
            return ZZ(v)

        def compute():
            ratio = ZZ(self.K.defining_polynomial().discriminant()) / self.discriminant
            return ZZ(ratio.abs().sqrt())
        return _try(compute)

    @cached_property
    def num_subfields(self):
        """Proper subfields, excluding Q and K itself."""
        v = self.known.get("num_subfields")
        if v is not None:
            return ZZ(v)
        return _try(lambda: ZZ(len(self.K.subfields()) - 2))

    @cached_property
    def is_galois(self):
        v = self.known.get("is_galois")
        if v is not None:
            return ZZ(1 if v else 0)
        return _try(lambda: ZZ(1 if self.K.is_galois() else 0))

    @cached_property
    def galois_label(self):
        """LMFDB label ``nTt`` of the Galois group, via PARI's ``polgalois``."""
        v = self.known.get("galois_label")
        if v is not None:
            return str(v)

        def compute():
            from sage.all import pari
            data = pari(self.K.defining_polynomial()).polgalois()
            return f"{self.degree}T{int(data[2])}"
        return _try(compute)

    # -- coordinates --------------------------------------------------------

    @cached_property
    def _basis_inv(self):
        return Matrix(QQ, [self.K(b).vector() for b in self.basis]).inverse()

    def coordinates(self, x):
        """Coefficient vector of ``x`` with respect to :attr:`basis`."""
        c = self.K(x).vector() * self._basis_inv
        if not all(t in ZZ for t in c):
            raise ValueError(f"{x} is not integral for the chosen basis")
        return tuple(ZZ(t) for t in c)

    def from_coordinates(self, coords):
        return sum(ZZ(c) * b for c, b in zip(coords, self.basis))

    # -- exact signs --------------------------------------------------------

    def signature(self, x):
        """
        Sign vector of ``x`` under all real embeddings, computed exactly in AA.

        Replaces the old ``NumberFieldData._signature``, which called
        ``.embeddings`` on an *element* (no such method) and then ``.sign()``,
        which needs the field to carry an embedding.
        """
        x = self.K(x)
        return tuple(ZZ(1) if e(x) > 0 else (ZZ(0) if e(x) == 0 else ZZ(-1))
                     for e in self.embeddings_exact)

    def is_totally_positive(self, x):
        return self.K(x).is_totally_positive()

    # -- units --------------------------------------------------------------

    @cached_property
    def fundamental_units(self):
        return [self.K(u) for u in self.K.unit_group().fundamental_units()]

    @cached_property
    def totally_positive_unit_basis(self):
        from .normalize import totally_positive_unit_basis
        return totally_positive_unit_basis(self.K)

    @cached_property
    def unit_signature_rank(self):
        """Rank over F_2 of the signature map on the unit group."""
        from sage.all import GF
        rows = [[ZZ(0) if s > 0 else ZZ(1) for s in self.signature(u)]
                for u in [self.K(-1)] + self.fundamental_units]
        return ZZ(Matrix(GF(2), rows).rank())

    @cached_property
    def normalizer(self):
        """The canonical-representative machinery, sharing this basis."""
        from .normalize import UnitOrbitNormalizer
        return UnitOrbitNormalizer(self.K, basis=self.basis,
                                   unit_basis=self.totally_positive_unit_basis,
                                   prec=self.prec)

    # -- geometry (inexact, candidate generation only) ----------------------

    def real_embeddings(self, prec=None):
        return self.K.real_embeddings(prec or self.prec)

    @cached_property
    def log_unit_lattice(self):
        """
        LLL-reduced basis of the log lattice of U^+, at the working precision.

        Reduced because the covering radius bound below is only useful for a
        reduced basis, and because it sits in an exponent in the search bound.
        LLL needs integers, so a scaled rational approximation is reduced and
        the resulting transformation applied to the exact real matrix.
        """
        R = RealField(self.prec)
        embs = self.real_embeddings()
        B = Matrix(R, [[R(e(u)).log() for e in embs]
                       for u in self.totally_positive_unit_basis])
        if B.nrows() == 0:
            return B
        scale = ZZ(2) ** (self.prec // 2)
        Bz = Matrix(ZZ, [[(x * scale).round() for x in row] for row in B.rows()])
        _, T = Bz.LLL(transformation=True)
        return T.change_ring(R) * B

    @cached_property
    def log_unit_covering_radius(self):
        """
        Upper bound on the covering radius of the log lattice of U^+.

        Babai's nearest-plane algorithm leaves an error ``e`` with
        ``|<e, b_i*>| <= ||b_i*||^2 / 2`` for each Gram-Schmidt vector, so
        ``rho <= sqrt(sum_i ||b_i*||^2) / 2``.

        Sage's ``gram_schmidt`` only supports ``RDF``/``CDF``, not
        ``RealField(prec)``, so the orthogonalisation is done here.  Dropping to
        ``RDF`` would work but would cap this at 53 bits, and it appears in an
        *exponent* in :attr:`trace_bound` -- a small error here is amplified.
        """
        R = RealField(self.prec)
        B = self.log_unit_lattice
        if B.nrows() == 0:
            return R(0)

        total, ortho = R(0), []
        for v in B.rows():
            w = v
            for u, sq in ortho:
                w = w - (v.dot_product(u) / sq) * u
            sq = w.dot_product(w)
            if sq <= 0:
                raise ValueError(
                    f"{self}: the log lattice of U^+ is rank deficient at "
                    f"{self.prec} bits; raise prec")
            ortho.append((w, sq))
            total += sq
        return total.sqrt() / 2

    @cached_property
    def trace_bound(self):
        """
        Rigorous cap on the trace of a canonical indecomposable.

        Uses N(alpha) <= |disc| (Kala-Yatsyna) together with the balancedness of
        the minimal-trace representative.  Strictly better than the classical
        regulator-dependent bound.
        """
        R = RealField(self.prec)
        n = self.degree
        return R(n) * R(self.discriminant.abs()) ** (R(1) / n) \
            * self.log_unit_covering_radius.exp()

    def __repr__(self):
        return f"<FieldContext {self.label or self.K} deg={self.degree}>"
