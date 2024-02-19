#!/usr/bin/env python

import os
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
