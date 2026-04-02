.PHONY: help install dev test test-minimal test-sage lint format clean

help:
	@echo "Indecomposables Project - Development Commands"
	@echo ""
	@echo "Setup:"
	@echo "  make install      Install core dependencies"
	@echo "  make dev          Install development dependencies"
	@echo ""
	@echo "Testing:"
	@echo "  make test         Run all tests (requires Sage)"
	@echo "  make test-minimal Run fast tests (no Sage required)"
	@echo "  make test-sage    Run full test suite with Sage"
	@echo "  make coverage     Generate coverage report"
	@echo ""
	@echo "Code Quality:"
	@echo "  make lint         Run flake8 linter"
	@echo "  make format       Auto-format code with black"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean        Remove cache and built files"

install:
	@echo "Install via conda/apt for Sage support"
	@echo "Ubuntu/Debian: sudo apt install sagemath"
	@echo "Or: conda install -c conda-forge sage"

dev:
	pip install -r requirements-dev.txt

test:
	python -m pytest tests/ -v

test-minimal:
	python -m pytest tests/test_main.py::TestMinimal -v

test-sage:
	python -m pytest tests/ -v -m "not slow"

coverage:
	python -m pytest tests/ --cov=scripts --cov-report=html --cov-report=term

lint:
	flake8 scripts/ tests/ --max-line-length=120 --count

format:
	black scripts/ tests/ --line-length=120

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".coverage" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "htmlcov" -exec rm -rf {} + 2>/dev/null || true
	@echo "Cleaned up cache files"

.DEFAULT_GOAL := help
