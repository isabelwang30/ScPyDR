import os
from setuptools import setup, find_packages

# Version-keeping code
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'scPyDR/version.py'), 'w') as fout:
    fout.write(
        "\n".join(["",
                   "# THIS FILE IS GENERATED FROM SETUP.PY",
                   "version = '{version}'",
                   "__version__ = version"]).format(version=VERSION)
    )

# Setup function
setup(
    name='scPyDR',
    version=VERSION,
    description='CSE185 Final Project',
    author='Anushka Sheoran, Isabel Wang, Monica Park',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "scpydr=scPyDR.scPyDR.__main__:main"  # Entry point for command-line script
        ],
    },
)
