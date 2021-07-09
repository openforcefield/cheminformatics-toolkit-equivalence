
from setuptools import setup
import pathlib

LONG_DESCRIPTION = pathlib.Path("README.md").read_text()

# This installer exists to make it easier to use off_coverage programmatically.

setup(
    name = "off_coverage",
    version = "1.0",
    description = "coverage reduction tool for the Open Force Field toolkit",
    
    long_description = LONG_DESCRIPTION,
    long_description_content_type = "text/markdown",
    
    author = "Andrew Dalke",
    author_email = 'dalke@dalkescientific.com',
    url = "https://github.com/openforcefield/cheminformatics-toolkit-equivalence/off_coverage",
    license = "MIT",
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
        ],

    py_modules = ["off_coverage"],

    ## Can't install as an entry point because this is a module, not a package.
    ## You'll need to use 'python -m off_coverage'.
    
    ## entry_points = {
    ##     "console_scripts": [
    ##         "off_coverage = off_coverage.main",
    ##         ],
    ##     },
    )
