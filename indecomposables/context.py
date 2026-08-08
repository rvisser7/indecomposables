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

ZZx = PolynomialRing(ZZ, "x")

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
        return float(v) if v is not None else float(self.K.regulator())

    @cached_property
    def class_number(self):
        v = self.known.get("class_number")
        return ZZ(v) if v is not None else ZZ(self.K.class_number())

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
    def log_unit_covering_radius(self):
        """
        Upper bound on the covering radius of the log lattice of U^+.

        Feeds the trace bound  Tr(alpha) <= n * disc^(1/n) * e^rho  for the
        minimal-trace representative of an indecomposable.
        """
        R = RealField(self.prec)
        embs = self.real_embeddings()
        B = Matrix(R, [[R(e(u)).log() for e in embs]
                       for u in self.totally_positive_unit_basis])
        if B.nrows() == 0:
            return R(0)
        gs, _ = B.gram_schmidt()
        return sum((row.norm() for row in gs.rows()), R(0)) / 2

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
