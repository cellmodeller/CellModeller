#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(name='CellModeller',
      packages=['CellModeller',
                'CellModeller.Biophysics',
                'CellModeller.Biophysics.BacterialModels',
                'CellModeller.Integration',
                'CellModeller.Regulation',
                'CellModeller.Signalling',
                'CellModeller.GUI'],
      package_data={'':['*.cl','*.ui']},
)
