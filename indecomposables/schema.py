"""
The database schema, defined once.

Everything downstream derives from ``COLUMNS``: the pipe-delimited row format,
the psycodict header/type preamble, the parser, ``data/README.md``, and the
split between the two files.  Nothing else in the codebase should hardcode a
column name, a column order, or a Postgres type.

Two files per degree, **one row per field in both**:

* ``data/degree<n>_mini.txt`` -- scalars only, complete coverage
* ``data/degree<n>.txt``      -- the same scalars plus every array column

``degree<n>_mini.txt`` is generated from ``degree<n>.txt`` by column projection
in ``scripts/merge_shards.py``, never computed separately, so the shared columns
cannot drift.

Per-signature results are stored as **arrays indexed by the ``signatures``
column**, not as extra rows -- which is what keeps mini a clean projection.
``signatures`` is ordered by :mod:`indecomposables.signatures`, so index 0 is
always the totally positive class.

Array types
-----------
Postgres array literals must be **rectangular**: ``{{1,2},{3,4,5}}` is not a
valid ``numeric[]``.  Columns whose nesting is genuinely jagged -- a different
number of minimal elements in each signature class, a different number of
vertices on each facet -- are therefore ``jsonb``, encoded as JSON text.
Rectangular columns stay as native arrays so they index and query normally.

Elements are stored as **coefficient vectors with respect to the integral basis
named in the ``basis`` column**, never as printed Sage expressions: printed forms
depend on the variable name and the defining polynomial, need a Sage ``eval`` to
read back, and are about twice the size.

Where the schema lives
----------------------
The column *declarations* are data and live in ``schema.yaml``; the codecs are
behaviour and live here, keyed by a small closed vocabulary of names.  Adding a
column is then a data-only diff that a non-Python consumer -- a Magma script, an
ingestion step, a front end -- can read directly, and a typo cannot break an
import.  What that costs is the safety net of Python refusing to load a
malformed file, which :func:`_validate` replaces with explicit checks at load
time.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field as _field
from pathlib import Path
from typing import Any, Callable

import yaml

__all__ = [
    "Column", "COLUMNS", "SCALAR_CODECS", "codec_for", "reserved", "load", "MINI", "FULL", "STATUSES",
    "columns", "header_lines", "encode_row", "decode_row", "project",
    "docs_table", "create_table_sql", "NULL", "SEP", "SchemaError",
]

SEP = "|"
NULL = r"\N"

MINI = "mini"
FULL = "full"

#: Every field in scope gets a row; ``status`` says what happened to it.
STATUSES = ("ok", "partial", "timeout", "error", "killed")


class SchemaError(ValueError):
    pass


# ---------------------------------------------------------------------------
# Codecs
# ---------------------------------------------------------------------------

def _enc_text(v):
    s = str(v)
    if SEP in s or "\n" in s:
        raise SchemaError(f"text value contains a delimiter: {s!r}")
    return s


def _enc_int(v):
    return str(int(v))


def _enc_float(v):
    return repr(float(v))


def _nested_encode(v, depth):
    """
    Encode a nested integer array.

    ``None`` inside an array becomes the literal ``NULL``.  Note this is *not*
    the same marker as a null whole field: ``\\N`` is COPY's field-level marker,
    while inside an array literal Postgres expects the word ``NULL``.

    Per-signature arrays need this: a class with no non-unit s-indecomposable
    has no minimum norm, and the entry has to say so rather than invent a
    sentinel that a reader could mistake for a real value.
    """
    if depth == 0:
        return "NULL" if v is None else str(int(v))
    return "{" + ",".join(_nested_encode(t, depth - 1) for t in v) + "}"


def _split_braces(s):
    if not (s.startswith("{") and s.endswith("}")):
        raise SchemaError(f"not an array literal: {s!r}")
    body, out, depth, cur = s[1:-1], [], 0, []
    if not body:
        return []
    for ch in body:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def _nested_decode(s, depth):
    if depth == 0:
        return None if s.upper() == "NULL" else int(s)
    return tuple(_nested_decode(t, depth - 1) for t in _split_braces(s))


def _array(depth):
    """Codec for a rectangular Postgres array of integers, nested ``depth`` deep."""
    return (lambda v: _nested_encode(v, depth), lambda s: _nested_decode(s, depth))


def _to_tuples(o):
    return tuple(_to_tuples(t) for t in o) if isinstance(o, (list, tuple)) else o


JSON = (lambda v: _enc_text(json.dumps(v, separators=(",", ":"))),
        lambda s: _to_tuples(json.loads(s)))

TEXT = (_enc_text, str)
INT = (_enc_int, int)
FLOAT = (_enc_float, float)
ARRAY1 = _array(1)
ARRAY2 = _array(2)
TEXT_ARRAY = (lambda v: "{" + ",".join(_enc_text(t) for t in v) + "}",
              lambda s: tuple(_split_braces(s)))


# ---------------------------------------------------------------------------
# Columns
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Column:
    name: str
    pg_type: str
    files: frozenset               # subset of {MINI, FULL}; empty means reserved
    doc: str                       # short_doc: the one-line description
    codec: tuple = _field(repr=False)
    nullable: bool = True
    long_doc: str = None           # optional prose for data/README.md

    @property
    def emitted(self):
        return bool(self.files)

    @property
    def encode(self) -> Callable[[Any], str]:
        return self.codec[0]

    @property
    def decode(self) -> Callable[[str], Any]:
        return self.codec[1]



#: Postgres type -> codec.  Every integer type shares one codec: Postgres stores
#: 1/0 for a smallint flag exactly as it stores 1/-1/0 for a tri-state, so a
#: separate boolean codec would only add an ambiguity the type cannot resolve.
SCALAR_CODECS = {
    "text": TEXT,
    "smallint": INT,
    "integer": INT,
    "bigint": INT,
    "numeric": INT,
    "double precision": FLOAT,
    "real": FLOAT,
    "jsonb": JSON,
}

SCHEMA_FILE = Path(__file__).with_name("schema.yaml")

#: Columns every row must carry.  Short enough to live here rather than in the
#: declarations, and it is really a property of the pipeline rather than of the
#: schema: without a label a row cannot be merged, sorted or joined.
NOT_NULL = ("lmfdb_label",)


def codec_for(pg_type):
    """
    Derive the codec from the declared Postgres type.

    Array dimension has to be written into the type -- Postgres calls both a 1-D
    and a 2-D array ``numeric[]``, so the type alone would not say how deeply to
    nest.  It accepts a multi-dimensional declaration and ignores the extra
    dimensions, so ``numeric[][]`` is valid DDL and records what we need.
    """
    base, dims = pg_type, 0
    while base.endswith("[]"):
        base, dims = base[:-2].strip(), dims + 1
    if dims == 0:
        if base not in SCALAR_CODECS:
            raise SchemaError(
                f"unknown type {pg_type!r}; known types are "
                f"{', '.join(sorted(SCALAR_CODECS))}, optionally with [] or [][]")
        return SCALAR_CODECS[base]
    if base not in SCALAR_CODECS or SCALAR_CODECS[base] is not INT:
        raise SchemaError(f"{pg_type!r}: only integer arrays are supported; "
                          "use jsonb for anything else")
    if dims == 1:
        return ARRAY1
    if dims == 2:
        return ARRAY2
    raise SchemaError(f"{pg_type!r}: arrays deeper than two dimensions must be "
                      "jsonb, since jagged nesting cannot be a Postgres array")


def _files(name, spec):
    """
    Which data files a column is written to.

    Required, and validated, because getting it wrong silently changes what is
    published.  ``mini`` without ``full`` is rejected: ``degree<n>_mini.txt`` is
    generated from ``degree<n>.txt`` by projection, so a mini-only column would
    have nothing to project from.  An empty list means the column is declared
    but not emitted -- a name and definition settled in advance of the data.
    """
    where = f"{SCHEMA_FILE.name}: column {name!r}"
    if "files" not in spec:
        raise SchemaError(
            f"{where}: missing required key 'files'; use [mini, full], [full], "
            "or [] for a column that is declared but written to no file")
    value = spec["files"] or []
    if isinstance(value, str):
        value = [value]
    bad = [v for v in value if v not in (MINI, FULL)]
    if bad:
        raise SchemaError(f"{where}: files has unknown name(s) {bad}; "
                          f"allowed values are {MINI!r} and {FULL!r}")
    if MINI in value and FULL not in value:
        raise SchemaError(
            f"{where}: files has {MINI!r} without {FULL!r}; the mini file is "
            "a projection of the full one, so a mini-only column cannot exist")
    return frozenset(value)


def _validate(name, spec):
    where = f"{SCHEMA_FILE.name}: column {name!r}"
    if not isinstance(spec, dict):
        raise SchemaError(f"{where}: expected a mapping of keys, got {type(spec).__name__}")
    for key in ("type", "short_doc"):
        if not spec.get(key):
            raise SchemaError(f"{where}: missing required key {key!r}")
    unknown = set(spec) - {"type", "files", "short_doc", "long_doc"}
    if unknown:
        raise SchemaError(
            f"{where}: unknown key(s) {sorted(unknown)}; "
            "allowed keys are type, tier, short_doc, long_doc")


def load(path=None):
    """
    Read the declarations.  Called once at import; re-callable for tests.

    The file is a flat mapping from column name to declaration, and its order is
    the column order in the data files -- YAML mappings load in document order,
    and :func:`header_lines` depends on that.
    """
    path = Path(path) if path else SCHEMA_FILE
    with path.open() as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict):
        raise SchemaError(f"{path}: expected a mapping from column name to declaration")

    cols = []
    for name, spec in doc.items():
        _validate(name, spec)
        pg = str(spec["type"]).strip()
        cols.append(Column(
            name=name, pg_type=pg,
            files=_files(name, spec),
            doc=" ".join(str(spec["short_doc"]).split()),
            long_doc=" ".join(str(spec.get("long_doc", "")).split()) or None,
            codec=codec_for(pg),
            nullable=name not in NOT_NULL))
    if not cols:
        raise SchemaError(f"{path}: no columns declared")
    return tuple(cols)


COLUMNS = load()

_BY_NAME = {c.name: c for c in COLUMNS}
if len(_BY_NAME) != len(COLUMNS):
    raise SchemaError("duplicate column name in COLUMNS")


def columns(tier: str) -> tuple:
    """
    Columns written to one file, in file order.

    Columns whose ``files`` is empty are declared but written nowhere, and so
    appear in neither.  ``full`` is a superset of ``mini`` by construction.
    """
    if tier not in (MINI, FULL):
        raise SchemaError(f"unknown file {tier!r}")
    return tuple(c for c in COLUMNS if tier in c.files)


def reserved() -> tuple:
    """Columns declared in the schema but not written to any file."""
    return tuple(c for c in COLUMNS if not c.emitted)


# ---------------------------------------------------------------------------
# Rows and files
# ---------------------------------------------------------------------------

def header_lines(tier: str) -> list:
    """psycodict preamble: names, Postgres types, blank line."""
    cols = columns(tier)
    return [SEP.join(c.name for c in cols),
            SEP.join(c.pg_type for c in cols),
            ""]


def encode_row(record, tier: str) -> str:
    out = []
    for c in columns(tier):
        v = record.get(c.name)
        if v is None:
            if not c.nullable:
                raise SchemaError(f"column {c.name!r} may not be NULL")
            out.append(NULL)
        else:
            out.append(c.encode(v))
    return SEP.join(out)


def decode_row(line: str, tier: str) -> dict:
    cols = columns(tier)
    parts = line.rstrip("\n\r").split(SEP)
    if len(parts) != len(cols):
        raise SchemaError(f"{tier} row has {len(parts)} fields, expected {len(cols)}")
    rec = {}
    for c, raw in zip(cols, parts):
        if raw == NULL:
            if not c.nullable:
                raise SchemaError(f"column {c.name!r} may not be NULL")
            rec[c.name] = None
        else:
            rec[c.name] = c.decode(raw)
    return rec


def project(line: str) -> str:
    """
    Turn one ``full`` line into the matching ``mini`` line.

    This is how ``degree<n>_mini.txt`` is produced: by projection, never by a
    second computation, so the shared columns cannot drift apart.
    """
    full = columns(FULL)
    parts = line.rstrip("\n\r").split(SEP)
    if len(parts) != len(full):
        raise SchemaError(f"full row has {len(parts)} fields, expected {len(full)}")
    keep = {c.name for c in columns(MINI)}
    return SEP.join(p for c, p in zip(full, parts) if c.name in keep)


def docs_table(tier: str) -> str:
    """Markdown for ``data/README.md``.  Generated; never hand-edit."""
    rows = ["| Column | Type | Description |", "| --- | --- | --- |"]
    for c in columns(tier):
        rows.append(f"| {c.name} | {c.pg_type} | {c.doc} |")
    return "\n".join(rows)


def create_table_sql(table: str, tier: str) -> str:
    """``CREATE TABLE`` for psycodict."""
    body = ",\n".join(f"    {c.name} {c.pg_type}" for c in columns(tier))
    return f"CREATE TABLE {table} (\n{body}\n);"


if __name__ == "__main__":
    for tier in (MINI, FULL):
        print(f"## `degree<n>{'_mini' if tier == MINI else ''}.txt`\n")
        print(docs_table(tier))
        print()
