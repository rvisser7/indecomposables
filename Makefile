# Short forms of the commands CI runs.  Everything here is also a plain
# command you can type; this is a convenience, not a build system.
#
# tox is deliberately absent.  Its selling points -- isolated virtualenvs and a
# matrix over Python versions -- do not apply when Sage pins the interpreter and
# cannot be pip-installed, which is why LMFDB's own tox.ini sets skipsdist=True
# and shells out to `sage`.  Test configuration lives in pyproject.toml instead.

SAGE ?= sage

# `pytest` rather than `python -m pytest`: conftest.py at the repo root puts
# the package on sys.path, so both work identically in a fresh clone.
PYTEST_FAST  = pytest tests/ -q --timeout=120
PYTEST_SAGE  = $(SAGE) -python -m pytest tests/ -q -n auto --timeout=600

.PHONY: help test test-sage test-doc test-slow verify lint docs validate all

help:
	@grep -E '^[a-z-]+:.*?## ' $(MAKEFILE_LIST) | \
	 awk -F':.*?## ' '{printf "  %-12s %s\n", $$1, $$2}'

test: ## Sage-free tests (seconds)
	$(PYTEST_FAST)

test-sage: ## full test suite under Sage
	$(PYTEST_SAGE)

test-doc: ## run BOTH doctest styles; each runner ignores the other's prompts
	$(SAGE) -t indecomposables/
	$(SAGE) -python -m pytest --doctest-modules indecomposables -q

test-slow: ## the slow tier (minutes)
	$(SAGE) -python -m pytest tests/ -q -m slow -n auto --durations=0

verify: ## cross-validate families against brute force (hours)
	$(SAGE) -python -m pytest tests/ -q -m verification -n auto --durations=0

lint: ## flake8, the same selection CI uses
	flake8 --count --select=E9,F63,F7,F82,F401,F841 --max-line-length=100 \
	       indecomposables scripts tests run_parallel.py

docs: ## regenerate data/README.md and data/schema.sql
	python scripts/make_docs.py

validate: ## check the committed data files
	python scripts/validate_data.py --all

all: lint test docs validate ## everything that does not need Sage
