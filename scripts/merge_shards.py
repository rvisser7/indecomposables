#!/usr/bin/env python3
"""
Merge worker shards into the published data files.

    python scripts/merge_shards.py --degree 3

Reads ``work/degree<n>/shard-*.txt`` (rows for fields that finished) together
with ``failures.txt`` and ``hard-timeouts.txt`` (labels that did not), and
produces:

    data/degree<n>.txt            sorted, deduplicated, with the psycodict preamble
    data/degree<n>_mini.txt       generated from the above by column projection
    data/degree<n>_failures.txt   label|reason for anything that did not finish
    data/manifest.json            scope, counts, schema version, sha256

Every column in the data files is an invariant of the field itself.  Facts about
the run -- code version, wall time, precision, the basis and order used -- go in
``manifest.json``, once per file rather than repeated on every row.  Fields that
did not finish are not rows at all, so no column needs a "not computed" NULL and
a gap is still accounted for: rows plus failures equals the scope.

``degree<n>_mini.txt`` is **never computed separately**.  It is a projection of
the full file, so the columns the two share cannot drift apart.

A label appearing in two shards is just a recomputation -- the later row wins.

No Sage is needed here; merging is pure text handling.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from indecomposables import schema                                  # noqa: E402

LABEL_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)\.(\d+)$")


def label_key(label: str):
    """Sort key: by absolute discriminant, then LMFDB index within it."""
    m = LABEL_RE.match(label)
    if not m:
        return (10 ** 18, 0, label)
    _, _, disc, index = m.groups()
    return (int(disc), int(index), label)


def read_shards(work_dir: Path):
    """All full-tier rows from every shard, as ``{label: record}``."""
    best, seen, malformed = {}, 0, []
    for shard in sorted(work_dir.glob("shard-*.txt")):
        for lineno, line in enumerate(shard.read_text().splitlines(), 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            seen += 1
            try:
                rec = schema.decode_row(line, schema.FULL)
            except schema.SchemaError as exc:
                malformed.append((shard.name, lineno, str(exc)))
                continue
            # Every shard row is a completed field, so a duplicate is just a
            # recomputation; the later one wins.
            best[rec["lmfdb_label"]] = rec
    return best, seen, malformed


def read_failures(work_dir: Path):
    """
    Fields that did not finish: ``{label: reason}``.

    Collected from the workers' own ``failures.txt`` and from the parent's
    ``hard-timeouts.txt`` (a watchdog-killed worker cannot write its own line).
    These become ``data/degree<n>_failures.txt`` rather than rows, so every row
    in the data file is a completed field and no column needs a "not computed"
    NULL.
    """
    out = {}
    for name in ("failures.txt", "hard-timeouts.txt"):
        path = work_dir / name
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split("|")
            if LABEL_RE.match(parts[0]):
                out[parts[0]] = parts[-1] if len(parts) > 1 else "unspecified"
    return out


def write_table(path: Path, records, tier):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = schema.header_lines(tier)
    lines += [schema.encode_row(r, tier) for r in records]
    path.write_text("\n".join(lines) + "\n")
    return path


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_provenance():
    """
    The per-run facts that used to be per-row columns.

    Recorded once per file.  A merge produces one file from one run, so this
    loses nothing unless rows from different code versions are merged into a
    single file -- which regenerating avoids.
    """
    def _run(cmd):
        try:
            import subprocess
            return subprocess.check_output(cmd, stderr=subprocess.DEVNULL,
                                           text=True).strip()
        except Exception:
            return None
    prov = {"code_version": _run(["git", "rev-parse", "--short", "HEAD"]),
            "generated": datetime.datetime.now(datetime.timezone.utc)
                                  .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "basis": "lmfdb_zk", "order": "maximal"}
    try:
        from sage.env import SAGE_VERSION
        prov["sage_version"] = SAGE_VERSION
    except Exception:
        pass
    return prov


def update_manifest(manifest_path: Path, entries):
    data = {}
    if manifest_path.exists():
        data = json.loads(manifest_path.read_text())
    for name, entry in entries.items():
        data[name] = entry
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--degree", type=int, required=True)
    p.add_argument("--work-dir", default="work")
    p.add_argument("--data-dir", default="data")
    p.add_argument("--selection-rule", default=None,
                   help="JSON string recording why this scope; stored in the manifest")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    work_dir = Path(args.work_dir) / f"degree{args.degree}"
    if not work_dir.is_dir():
        raise SystemExit(f"no shards at {work_dir}")

    records, seen, malformed = read_shards(work_dir)
    superseded = seen - len(records)
    failures = {k: v for k, v in read_failures(work_dir).items() if k not in records}

    ordered = [records[k] for k in sorted(records, key=label_key)]

    print(f"degree {args.degree}")
    print(f"  shard rows read : {seen}")
    print(f"  unique labels   : {len(ordered)}  ({superseded} superseded by a better run)")
    print(f"  did not finish  : {len(failures)}")
    if malformed:
        print(f"  MALFORMED       : {len(malformed)}")
        for name, lineno, msg in malformed[:5]:
            print(f"    {name}:{lineno}  {msg}")
    if args.dry_run:
        return 0

    data_dir = Path(args.data_dir)
    full_path = write_table(data_dir / f"degree{args.degree}.txt", ordered, schema.FULL)

    # mini by projection, so the shared columns cannot disagree
    full_lines = full_path.read_text().splitlines()
    mini_lines = schema.header_lines(schema.MINI) + \
        [schema.project(line) for line in full_lines[3:] if line]
    mini_path = data_dir / f"degree{args.degree}_mini.txt"
    mini_path.write_text("\n".join(mini_lines) + "\n")

    discs = [r["discriminant"] for r in ordered if r.get("discriminant") is not None]
    if failures:
        fail_path = data_dir / f"degree{args.degree}_failures.txt"
        fail_path.write_text("lmfdb_label|reason\n" + "".join(
            f"{k}|{failures[k]}\n" for k in sorted(failures, key=label_key)))
        print(f"  wrote {fail_path} ({len(failures)} field(s))")

    entry = dict(run_provenance())
    entry.update({"degree": args.degree, "rows": len(ordered),
                  "failures": len(failures),
                  "disc_min": min(discs) if discs else None,
                  "disc_max": max(discs) if discs else None})
    if args.selection_rule:
        entry["selection_rule"] = json.loads(args.selection_rule)
    update_manifest(data_dir / "manifest.json", {
        full_path.name: dict(entry, tier="full", sha256=sha256(full_path)),
        mini_path.name: dict(entry, tier="mini", sha256=sha256(mini_path)),
    })

    print(f"  wrote {full_path} ({full_path.stat().st_size:,} bytes)")
    print(f"  wrote {mini_path} ({mini_path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
