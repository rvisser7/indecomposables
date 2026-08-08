"""
Structural checks on the package's own imports.

These run without Sage, so they catch a whole class of breakage in the fast CI
tier rather than at the first real computation.
"""

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parent.parent / "indecomposables"
MODULES = {p.stem for p in PACKAGE.glob("*.py")} | {"families"}


def _python_files():
    return sorted(PACKAGE.rglob("*.py"))


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: p.name)
def test_sibling_imports_are_relative(path):
    """
    A module inside the package must import its siblings relatively.

    ``from normalize import ...`` works when the package directory happens to be
    on sys.path -- which it is when you run a file directly, or from a notebook
    in that directory -- and fails everywhere else.  The symptom is a
    ``ModuleNotFoundError`` from deep inside a worker, long after import time,
    for a module that plainly exists.
    """
    tree = ast.parse(path.read_text())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            if node.module.split(".")[0] in MODULES:
                bad.append(f"line {node.lineno}: from {node.module} import ...")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] in MODULES:
                    bad.append(f"line {node.lineno}: import {alias.name}")
    assert not bad, (
        f"{path.name} imports a sibling absolutely; use a relative import:\n  "
        + "\n  ".join(bad))


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: p.name)
def test_no_stdlib_shadowing_names(path):
    """
    No module may be named after a standard library module.

    ``logging.py`` inside the package shadows the stdlib for any tool run from
    that directory, and the failure surfaces inside the tool rather than here.
    """
    import sys
    stdlib = getattr(sys, "stdlib_module_names", frozenset())
    assert path.stem not in stdlib, (
        f"{path.name} shadows the standard library module {path.stem!r}")
