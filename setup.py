#!/usr/bin/env python

import os
from setuptools import setup
import subprocess

# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
version_py = os.path.join(os.path.dirname(__file__), 'CellModeller', 'version.py')

def get_git_version(default="0.1.0"):
    try:
        # Ensure the output from git is decoded to a string.
        version_git = subprocess.check_output(["git", "describe"], text=True).strip()
        # Format the version string to be PEP 440 compliant.
        # Example: 'v4.3-45-g0b90179' -> '4.3.45+g0b90179'
        # Adjust the formatting as per your versioning scheme.
        version_git = version_git.lstrip('v').replace('-', '.').replace('.', '+git', 1)
    except Exception:
        # If git is not available, read the version from version.py
        with open(version_py, 'r') as fh:
            version_git_contents = fh.read().strip()
            version_git = version_git_contents.split('=')[-1].replace('"', '').strip()
    return version_git

# Use the function to get the version.
version = get_git_version()

setup(
    name='CellModeller',
    install_requires=['numpy', 'scipy', 'pyopengl', 'mako', 'pyqt5', 'pyopencl', 'reportlab', 'matplotlib'],
    setup_requires=['numpy', 'scipy', 'pyopengl', 'mako', 'pyqt5', 'pyopencl', 'reportlab', 'matplotlib'],
    packages=['CellModeller',
              'CellModeller.Biophysics',
              'CellModeller.Biophysics.BacterialModels',
              'CellModeller.Biophysics.GeneralModels',
              'CellModeller.Integration',
              'CellModeller.Regulation',
              'CellModeller.Signalling',
              'CellModeller.GUI'],
    package_data={'': ['*.cl', '*.ui']},
    python_requires='>=3',
    version=version
import re
from setuptools import setup
import subprocess

version_py = os.path.join(os.path.dirname(__file__), 'CellModeller/version.py')

try:
    version_git = subprocess.check_output(["git", "describe"], text=True).strip()
    # Transform git tag into a PEP 440 compliant version string
    version_git = re.sub(r'^v', '', version_git)
    version_git = re.sub(r'-(\d+)-g([0-9a-f]{7,})', r'+\1.sha\2', version_git)
except:
    with open(version_py, 'r') as fh:
        version_git = fh.read().strip().split('=')[-1].replace('"', '')

setup(
    name='CellModeller',
    install_requires=[
        'numpy', 'scipy', 'pyopengl', 'mako', 'pyqt5', 'pyopencl', 'reportlab', 'matplotlib'
    ],
    packages=[
        'CellModeller',
        'CellModeller.Biophysics',
        'CellModeller.Biophysics.BacterialModels',
        'CellModeller.Biophysics.GeneralModels',
        'CellModeller.Integration',
        'CellModeller.Regulation',
        'CellModeller.Signalling',
        'CellModeller.GUI'
    ],
    package_data={'': ['*.cl', '*.ui']},
    python_requires='>=3',
    version=str(version_git)
)

