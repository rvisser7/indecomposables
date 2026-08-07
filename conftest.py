"""
Root conftest.

Puts the repository root on ``sys.path`` so a bare ``pytest`` works in a fresh
clone, before anyone has run ``pip install -e .``.  Without this, ``pytest``
(the console script) cannot import the package, while ``python -m pytest``
can -- because that form silently prepends the working directory.  Two commands
that look equivalent behaving differently is exactly the kind of thing that
wastes a new contributor's first half hour.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
  
