"""
mm_project
A package for doing MC of a Lennard Jones fluid
"""
import os
import sys
from setuptools import setup, find_packages, Extension
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:]),


# Build our C++ module
# NOTE: Pybind11/Eigen were installed into CONDA_PREFIX
#       so we need to add that to the include paths
conda_prefix = os.environ['CONDA_PREFIX']
eigen_path = os.path.join(conda_prefix, 'include', 'eigen3')
qm_cpp_module = Extension('mm_project.mm_cpp',
                          include_dirs = [eigen_path],
                          sources = ['mm_cpp/day_6.cpp',
                                     'mm_cpp/export.cpp'])


setup(
    # Self-descriptive entries which should always be present
    name='mm_project',
    author='The Molecular Sciences Software Institute',
    author_email='janash@vt.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Include the compiled extension
    ext_modules = [qm_cpp_module]

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
