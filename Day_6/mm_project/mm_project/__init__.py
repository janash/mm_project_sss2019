"""
mm_project
A package for doing MC of a Lennard Jones fluid
"""

# Add imports here
from .monte_carlo import *
from .coordinates import generate_initial_coordinates, Box

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
