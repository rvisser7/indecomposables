#!/usr/bin/env python3
"""
Check the published data files.

    python scripts/validate_data.py --degree 3
    python scripts/validate_data.py --all

Two kinds of check.  **Structural** ones verify the file matches the schema:
preamble, arity, types, label ordering, uniqueness, and that the mini file is
exactly the projection of the full one.  **Mathematical** ones verify internal
consistency of quantities computed by different routes -- these are the ones
that catch real bugs.

The most valuable single check is

    len(signatures) == 2^(degree - unit_signature_rank) == narrow_class_number / class_number

which ties together three columns arrived at completely independently: a coset
enumeration, an F_2 rank, and a pair of LMFDB class numbers.  A wrong embedding
order, a botched U^+ computation, or a mis-joined LMFDB row all break it.

Exits non-zero if anything fails, so it can be a CI step.  Needs no Sage.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from indecomposables import schema                                  # noqa: E402
from scripts.merge_shards import label_key, LABEL_RE                # noqa: E402


class Report:
    def __init__(self):
        self.errors, self.warnings, self.checks = [], [], 0

    def check(self, ok, where, msg):
        self.checks += 1
        if not ok:
            self.errors.append(f"{where}: {msg}")

    def warn(self, ok, where, msg):
        self.checks += 1
        if not ok:
            self.warnings.append(f"{where}: {msg}")


def read_table(path: Path, tier):
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path}: too short to contain a preamble")
    cols = [c.name for c in schema.columns(tier)]
    types = [c.pg_type for c in schema.columns(tier)]
    if lines[0].split("|") != cols:
        raise ValueError(f"{path}: header does not match schema tier {tier!r}")
    if lines[1].split("|") != types:
        raise ValueError(f"{path}: type line does not match schema tier {tier!r}")
    if lines[2].strip():
        raise ValueError(f"{path}: third line must be blank")
    return [schema.decode_row(ln, tier) for ln in lines[3:] if ln.strip()], lines


def validate_degree(degree, data_dir: Path, rep: Report, required=False):
    full_path = data_dir / f"degree{degree}.txt"
    mini_path = data_dir / f"degree{degree}_mini.txt"
    if not full_path.exists():
        # An explicitly requested degree that is missing is an error: it usually
        # means a merge did not run.  A degree only picked up by --all is not.
        if required:
            rep.errors.append(f"{full_path}: missing")
        else:
            rep.warnings.append(f"{full_path}: missing, skipped")
        return

    rows, full_lines = read_table(full_path, schema.FULL)
    where = full_path.name

    # ---- structure -------------------------------------------------------
    labels = [r["lmfdb_label"] for r in rows]
    rep.check(len(labels) == len(set(labels)), where, "duplicate labels")
    rep.check(labels == sorted(labels, key=label_key), where,
              "rows are not sorted by (|disc|, index)")

    if mini_path.exists():
        mini_lines = mini_path.read_text().splitlines()
        expected = schema.header_lines(schema.MINI) + \
            [schema.project(ln) for ln in full_lines[3:] if ln.strip()]
        rep.check([ln for ln in mini_lines if ln.strip() or True] == expected, where,
                  "mini file is not the projection of the full file "
                  "(regenerate with merge_shards.py, do not hand-edit)")
    else:
        rep.warnings.append(f"{mini_path}: missing")

    # ---- per row ---------------------------------------------------------
    for r in rows:
        lab = r["lmfdb_label"]
        w = f"{where}[{lab}]"
        m = LABEL_RE.match(lab or "")
        rep.check(bool(m), w, f"malformed label {lab!r}")
        if m:
            rep.check(int(m.group(1)) == degree, w, "label degree disagrees with file")
            rep.check(int(m.group(1)) == int(m.group(2)), w,
                      "label signature is not (n,0), so the field is not totally real")
            rep.check(abs(r["discriminant"]) == int(m.group(3)), w,
                      "discriminant column disagrees with the label")
        if r.get("degree") is not None:      # full file only
            rep.check(r["degree"] == degree, w,
                      "degree column disagrees with file")
        # Every row is a completed field: unfinished ones are listed in
        # degree<n>_failures.txt rather than carried here as NULLs.
        rep.check(r.get("num_indecomposables_all") is not None, w,
                  "row has no results; unfinished fields belong in the failures file")

        _check_per_signature_lengths(r, w, rep)
        _check_signatures(r, w, rep)
        _check_indecomposables(r, w, rep)
        _check_sail(r, w, rep)

    print(f"  degree {degree}: {len(rows)} rows")


#: Columns holding one value per signature class.  Derived rather than listed,
#: so a new per-signature column is checked automatically.
def _per_signature(rec):
    fixed = {"coeffs", "fundamental_units", "positive_unit_basis", "signature_classes"}
    return [c for c in schema.columns(schema.FULL)
            if c.name not in fixed
            and (c.pg_type.endswith("[]") or c.pg_type == "jsonb")
            and rec.get(c.name) is not None]


def _check_per_signature_lengths(r, w, rep):
    """
    Every per-signature array must have exactly num_signature_classes entries.

    Now that all indecomposable and sail data is indexed by signature class,
    this single check covers fifteen columns at once, and it is the check that
    catches a misalignment between element data and sail data -- the failure
    mode that would otherwise make sail_facet_vertices point at the wrong
    class's elements.
    """
    nclass = r.get("num_signature_classes")
    if nclass is None:
        return
    for col in _per_signature(r):
        rep.check(len(r[col.name]) == nclass, w,
                  f"{col.name} has {len(r[col.name])} entries, "
                  f"num_signature_classes says {nclass}")


def _check_signatures(r, w, rep):
    n, sigs = r.get("degree"), r.get("signature_classes")
    if n is None and sigs:
        n = len(sigs[0])
    rank, nclass = r.get("unit_signature_rank"), r.get("num_signature_classes")
    if sigs is None or rank is None:
        return
    rep.check(len(sigs) == nclass, w,
              f"signature_classes has {len(sigs)} entries, "
              f"num_signature_classes says {nclass}")
    rep.check(nclass == 2 ** (n - rank), w,
              f"num_signature_classes {nclass} != 2^({n}-{rank})")
    h, hp = r.get("class_number"), r.get("narrow_class_number")
    if h and hp:
        rep.check(hp % h == 0 and hp // h == nclass, w,
                  f"h+/h = {hp}/{h} != num_signature_classes {nclass}")
    rep.check(all(len(s) == n for s in sigs), w, "a signature vector has the wrong length")
    rep.check(all(x in (1, -1) for s in sigs for x in s), w,
              "signature entries must be +-1")
    rep.check(sigs[0] == tuple([1] * n), w,
              "signature class 0 must be the totally positive one")
    masks = r.get("signature_masks")
    if masks is not None:
        rep.check(list(masks) == sorted(masks) and masks[0] == 0, w,
                  "signature_masks must be ascending and start at 0")
        rep.check(len(masks) == nclass, w, "signature_masks has the wrong length")
        for mask, sig in zip(masks, sigs):
            got = sum((1 << (n - 1 - i)) for i, s in enumerate(sig) if s < 0)
            rep.check(got == mask, w, f"mask {mask} does not encode {sig}")



def _check_indecomposables(r, w, rep):
    """
    Element-level checks.

    Shape checks apply to every signature class; the classical norm and trace
    bounds apply only to the totally positive class, index 0, which is where
    "indecomposable" means what Dress-Scharlau and Kala-Yatsyna mean by it.
    """
    ind = r.get("indecomposables_all")
    if ind is None:
        return
    degree = r.get("degree")

    # -- shape, every class ------------------------------------------------
    for name in ("indecomposables_norms_all", "indecomposables_traces_all",
                 "indecomposables_on_sail"):
        arr = r.get(name)
        if arr is not None:
            rep.check([len(x) for x in arr] == [len(x) for x in ind], w,
                      f"{name} lengths disagree with indecomposables_all")
    if r.get("num_indecomposables_all") is not None:
        rep.check([len(x) for x in ind] == list(r["num_indecomposables_all"]), w,
                  "num_indecomposables_all disagrees with indecomposables_all")
    if degree:
        rep.check(all(len(v) == degree for cls in ind for v in cls), w,
                  "a coefficient vector has the wrong length")

    norms_all = r.get("indecomposables_norms_all")
    sail_all = r.get("indecomposables_on_sail")
    if sail_all is not None:
        rep.check(all(f in (0, 1) for cls in sail_all for f in cls), w,
                  "on_sail flags must be 0 or 1")
        if r.get("num_indecomposables_on_sail") is not None:
            rep.check([sum(x) for x in sail_all] == list(r["num_indecomposables_on_sail"]), w,
                      "num_indecomposables_on_sail disagrees with indecomposables_on_sail")

    # -- per class extremes ------------------------------------------------
    if norms_all is not None:
        rep.check(all(n >= 1 for cls in norms_all for n in cls), w,
                  "an absolute norm is not a positive integer")
        for name, want in (("min_norm_indecomposable",
                            [min([n for n in cls if n != 1], default=None)
                             for cls in norms_all]),
                           ("max_norm_indecomposable",
                            [max(cls, default=None) for cls in norms_all])):
            if r.get(name) is not None:
                rep.check(list(r[name]) == want, w,
                          f"{name} disagrees with indecomposables_norms_all "
                          "(NULL is correct for a class with no such element)")

    # -- the totally positive class ----------------------------------------
    if not norms_all:
        return
    tp = list(norms_all[0])
    if not tp:
        return
    rep.check(1 in tp, w,
              "1 is always indecomposable but the totally positive class has "
              "no norm-1 element")
    disc = abs(r["discriminant"])
    rep.check(max(tp) <= disc, w,
              f"norm {max(tp)} exceeds |disc| {disc} (Kala-Yatsyna)")
    if degree == 2:
        # Dress-Scharlau give the sharper N(alpha) <= disc/4 in degree 2
        rep.check(4 * max(tp) <= disc, w,
                  f"norm {max(tp)} exceeds disc/4 = {disc / 4} (Dress-Scharlau)")
        if sail_all is not None:
            rep.check(all(sail_all[0]), w,
                      "in degree 2 every indecomposable lies on the sail")
    if degree and r.get("indecomposables_traces_all"):
        rep.check(all(t >= degree for t in r["indecomposables_traces_all"][0]), w,
                  "a trace is below n, impossible for a totally positive element")


def _check_sail(r, w, rep):
    verts = r.get("sail_facet_vertices")
    mins = r.get("indecomposables_all")
    nfac = r.get("num_sail_facets_by_signature")
    if verts is None:
        return
    if nfac is not None:
        rep.check([len(f) for f in verts] == list(nfac), w,
                  "sail_facet_vertices disagrees with num_sail_facets_by_signature")
    if mins is not None:
        for ci, facets in enumerate(verts):
            for facet in facets:
                rep.check(all(0 <= i < len(mins[ci]) for i in facet), w,
                          f"a facet in class {ci} indexes outside indecomposables_all")
                rep.check(len(set(facet)) == len(facet), w,
                          f"a facet in class {ci} repeats a vertex")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--degree", type=int, action="append")
    p.add_argument("--all", action="store_true")
    p.add_argument("--data-dir", default="data")
    p.add_argument("-q", "--quiet", action="store_true")
    args = p.parse_args(argv)

    data_dir = Path(args.data_dir)
    explicit = bool(args.degree)
    degrees = args.degree or []
    if args.all or not degrees:
        degrees = sorted(int(p.stem.replace("degree", ""))
                         for p in data_dir.glob("degree*.txt")
                         if p.stem.replace("degree", "").isdigit())
    if not degrees:
        if args.all:
            # A repository with no data yet is not a failure; only an explicit
            # --degree for a missing file is.
            print(f"no degree files in {data_dir}, nothing to validate")
            return 0
        raise SystemExit(f"no degree files found in {data_dir}")

    rep = Report()
    print(f"validating {data_dir}")
    for d in degrees:
        validate_degree(d, data_dir, rep, required=explicit)

    print(f"\n{rep.checks} checks run")
    for msg in rep.warnings:
        print(f"  WARN  {msg}")
    for msg in rep.errors[:50]:
        print(f"  FAIL  {msg}")
    if len(rep.errors) > 50:
        print(f"  ... and {len(rep.errors) - 50} more")
    if rep.errors:
        print(f"\n{len(rep.errors)} failure(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
