"""
Logger setup.

Deliberately *not* named ``logging.py``: a module of that name inside the
package directory shadows the standard library for any tool run from there
(flake8, pytest with rootdir on the path), which fails in confusing ways.

Two changes from the original:

* nothing is configured at import time.  ``main.py`` used to call
  ``setup_logger()`` at module scope, so importing the library reconfigured
  logging for whatever process imported it -- including pytest;
* the library only ever calls ``logging.getLogger(__name__)``.  Handlers are
  attached by the *application* (``run_parallel.py``, a script, a notebook),
  which is the one place that knows where output should go.
"""

from __future__ import annotations

import logging
import os

ROOT = "indecomposables"


def get_logger(name=None):
    """The logger a library module should use.  Attaches no handlers."""
    return logging.getLogger(ROOT if name is None else f"{ROOT}.{name}")


def setup_logging(log_file=None, verbose=False, debug=False):
    """
    Attach handlers.  Call once, from an application entry point.

    Console gets WARNING (or INFO with ``verbose``) as bare messages; the file,
    if given, gets INFO (or DEBUG) with timestamps.
    """
    logger = logging.getLogger(ROOT)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.WARNING)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    if log_file:
        parent = os.path.dirname(log_file)
        if parent:
            os.makedirs(parent, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG if debug else logging.INFO)
        fh.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(fh)
    return logger
