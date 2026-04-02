"""
Pytest configuration and shared fixtures for indecomposables tests.
"""

import pytest
import sys
import os

# Add scripts directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))


@pytest.fixture(scope="session")
def has_sage():
    """Check if Sage is available in environment."""
    try:
        import sage.all
        return True
    except ImportError:
        return False


@pytest.fixture(scope="session")
def sage_available(has_sage):
    """Skip test if Sage is not available."""
    if not has_sage:
        pytest.skip("Sage not available")


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "minimal: fast tests without Sage (for CI)"
    )
    config.addinivalue_line(
        "markers", "full: full tests with Sage"
    )
    config.addinivalue_line(
        "markers", "slow: slow computational tests"
    )
