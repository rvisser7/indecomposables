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

from dataclasses import dataclass, field as _field
from pathlib import Path
from typing import Any, Callable

import yaml

__all__ = [
    "Column", "COLUMNS", "SCALAR_CODECS", "codec_for", "pg_type_for", "parse_type", "reserved", "load", "MINI", "FULL", "STATUSES",
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


#: Delimiters for a nested value.  Square brackets, not Postgres's braces: the
#: data files are read on GitHub far more often than they are loaded, order
#: matters here, and braces read as set notation.  Converting to Postgres is a
#: mechanical substitution at load time.
OPEN, CLOSE = "[", "]"

#: Bracket pairs accepted when *reading*.  Only ``[]`` is ever written, but
#: earlier files used Postgres's ``{}``, and a format change should not force a
#: recomputation -- so the parser takes either.  Write strictly, read leniently.
BRACKETS = {"[": "]", "{": "}"}


def _enc_nested(v):
    """
    Encode a value of any nesting depth.

    One encoder for every nested column.  The separate rectangular-array and
    JSON codecs existed only because Postgres array literals must be rectangular
    and jagged data therefore had to be jsonb.  With square brackets the two are
    written identically, so the distinction survives only in the declared type,
    where it decides the Postgres column type rather than the file syntax.
    """
    if v is None:
        return NULL
    if isinstance(v, (list, tuple)):
        return OPEN + ",".join(_enc_nested(t) for t in v) + CLOSE
    return str(int(v))


def _split_list(s):
    """Split the body of a bracketed list on top-level commas."""
    opener = s[:1]
    if opener not in BRACKETS or not s.endswith(BRACKETS[opener]):
        raise SchemaError(f"not a list: {s!r}")
    body, out, depth, cur = s[1:-1], [], 0, []
    if not body:
        return []
    for ch in body:
        if ch in BRACKETS:
            depth += 1
        elif ch in BRACKETS.values():
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return out


def _dec_nested(s):
    """Inverse of :func:`_enc_nested`.  Depth is read off the brackets."""
    s = s.strip()
    if s == NULL:
        return None
    if s[:1] in BRACKETS:
        return tuple(_dec_nested(t) for t in _split_list(s))
    if s.upper() == "NULL":              # in-array null, as older files wrote it
        return None
    return int(s)


def _to_tuples(o):
    return tuple(_to_tuples(t) for t in o) if isinstance(o, (list, tuple)) else o


TEXT = (_enc_text, str)
INT = (_enc_int, int)
FLOAT = (_enc_float, float)
NESTED = (_enc_nested, _dec_nested)

#: Kept as names for readability; all nested columns share one codec now.
ARRAY1 = ARRAY2 = ARRAY3 = JSON = NESTED


@dataclass(frozen=True)
class Column:
    name: str
    declared_type: str             # logical: element type plus nesting depth
    jagged: bool                   # ragged nesting, so stored as jsonb
    pg_type: str                   # derived: what goes in the file preamble
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
}

SCHEMA_FILE = Path(__file__).with_name("schema.yaml")

#: Columns every row must carry.  Short enough to live here rather than in the
#: declarations, and it is really a property of the pipeline rather than of the
#: schema: without a label a row cannot be merged, sorted or joined.
NOT_NULL = ("lmfdb_label",)


def parse_type(declared):
    """
    Split a declared type into ``(element type, nesting depth)``.

    ``type`` is a *logical* description: the element type plus how deeply it is
    nested.  ``numeric[][]`` means a 2-D array of numbers whether or not
    Postgres can store it that way; :func:`pg_type_for` decides that.
    """
    base, dims = str(declared).strip(), 0
    while base.endswith("[]"):
        base, dims = base[:-2].strip(), dims + 1
    if base not in SCALAR_CODECS:
        raise SchemaError(
            f"unknown element type {base!r} in {declared!r}; known types are "
            f"{', '.join(sorted(SCALAR_CODECS))}")
    if dims > 3:
        raise SchemaError(f"{declared!r}: nesting deeper than three is not supported")
    if dims and SCALAR_CODECS[base] is not INT:
        raise SchemaError(f"{declared!r}: only integer element types can be nested")
    return base, dims


def pg_type_for(declared, jagged=False):
    """
    The type written into the data file preamble and ``CREATE TABLE``.

    A **jagged** nested column becomes ``jsonb``: Postgres array literals must be
    rectangular, so ``{{1,2},{3,4,5}}`` is rejected on COPY even though the
    ``CREATE TABLE`` would have been accepted.  A rectangular one keeps its array
    type -- Postgres accepts a multi-dimensional declaration and ignores the
    extra dimensions, so ``numeric[][]`` is valid DDL.
    """
    base, dims = parse_type(declared)
    if dims == 0:
        return base
    return "jsonb" if jagged else f"{base}{'[]' * dims}"


def codec_for(declared, jagged=False):
    """A scalar codec, or the single nested codec for anything with brackets."""
    base, dims = parse_type(declared)
    return SCALAR_CODECS[base] if dims == 0 else NESTED


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
    _, dims = parse_type(spec["type"])
    if spec.get("jagged") and dims == 0:
        raise SchemaError(f"{where}: 'jagged' only makes sense for a nested type")
    if dims == 3 and not spec.get("jagged"):
        raise SchemaError(
            f"{where}: a three-deep column must be declared jagged, since a "
            "rectangular three-dimensional Postgres array is not something we "
            "produce; add `jagged: true`")
    unknown = set(spec) - {"type", "files", "short_doc", "long_doc", "jagged"}
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
        declared = str(spec["type"]).strip()
        jagged = bool(spec.get("jagged", False))
        cols.append(Column(
            name=name, declared_type=declared, jagged=jagged,
            pg_type=pg_type_for(declared, jagged),
            files=_files(name, spec),
            doc=" ".join(str(spec["short_doc"]).split()),
            long_doc=" ".join(str(spec.get("long_doc", "")).split()) or None,
            codec=codec_for(declared, jagged),
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


def _shape(v, limit=3):
    """A short description of a value's nesting, for error messages."""
    depth, probe = 0, v
    while isinstance(probe, (list, tuple)) and depth < limit + 1:
        if not probe:
            return f"{'list of ' * depth}empty list"
        probe = probe[0]
        depth += 1
    return f"{'list of ' * depth}{type(probe).__name__}"


def _actual_depth(v):
    """Nesting depth of a value: 0 for a scalar, 1 for a flat list, and so on."""
    depth = 0
    while isinstance(v, (list, tuple)):
        depth += 1
        if not v:
            break                       # an empty list ends the descent
        v = v[0]
    return depth


def _check_shape(col, v):
    """
    Verify a value against its declared type before encoding.

    Two checks, both of which the old per-depth codecs gave for free and the
    single recursive encoder does not: it will happily write any nesting at all,
    so a producer emitting three levels into a two-level column would now
    succeed silently and corrupt the file.

    * **depth** must match the declared type, always;
    * **rectangularity** is required unless the column is declared ``jagged``,
      since that is what decides whether it can be a Postgres array on load.
    """
    _, dims = parse_type(col.declared_type)
    if dims == 0 or v is None:
        return

    if not v:
        return                          # an empty list is fine at any depth
    got = _actual_depth(v)
    if got != dims:
        raise SchemaError(
            f"column {col.name!r} is declared {col.declared_type!r} "
            f"(nesting {dims}) but the value has nesting {got}. The producer and "
            "schema.yaml disagree about this column's shape.")
    if col.jagged or dims < 2:
        return
    lengths = {len(row) for row in v if isinstance(row, (list, tuple))}
    if len(lengths) > 1:
        raise SchemaError(
            f"column {col.name!r} is declared {col.declared_type!r} but the data "
            f"is ragged (row lengths {sorted(lengths)}). Add `jagged: true` if "
            "that is intended -- it is what decides the Postgres column type.")


def encode_row(record, tier: str) -> str:
    out = []
    for c in columns(tier):
        v = record.get(c.name)
        if v is None:
            if not c.nullable:
                raise SchemaError(f"column {c.name!r} may not be NULL")
            out.append(NULL)
            continue
        try:
            _check_shape(c, v)
            out.append(c.encode(v))
        except SchemaError:
            raise
        except (TypeError, ValueError) as exc:
            # A bare TypeError from inside a codec gives no hint which column
            # is wrong.  Usually the producer and schema.yaml disagree about a
            # column's nesting -- for instance a per-signature column declared
            # as a two-deep array while the code emits three levels.
            raise SchemaError(
                f"column {c.name!r} declared as {c.pg_type!r} but got "
                f"{_shape(v)}: {exc}. If these disagree, the code and "
                "schema.yaml are out of sync -- check the column's `type`."
            ) from exc
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
            continue
        try:
            rec[c.name] = c.decode(raw)
        except SchemaError:
            raise
        except (TypeError, ValueError) as exc:
            # Without this the failure is a bare int() error from inside the
            # parser, with no indication of which column or which row.
            raise SchemaError(
                f"column {c.name!r} (declared {c.declared_type!r}) could not be "
                f"read from {raw[:60]!r}: {exc}") from exc
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
