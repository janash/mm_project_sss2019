"""
Unit and regression test for the mm_project package.
"""

# Import package, test suite, and other packages as needed
import mm_project
import pytest
import sys

def test_mm_project_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mm_project" in sys.modules
