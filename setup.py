#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(name='CellModeller',
      install_requires=['numpy', 'scipy', 'pyopengl', 'mako', 'pyopencl'],
      setup_requires=['numpy', 'scipy', 'pyopengl', 'mako', 'pyopencl'],
      packages=['CellModeller',
                'CellModeller.Biophysics',
                'CellModeller.Biophysics.BacterialModels',
                'CellModeller.Integration',
                'CellModeller.Regulation',
                'CellModeller.Signalling',
                'CellModeller.GUI'],
      package_data={'CellModeller.GUI':['PyGLViewer.ui'],
                    'CellModeller.Integration':['CLCrankNicIntegrator.cl'],
                    'CellModeller.Integration':['CLEulerIntegrator.cl'],
                    'CellModeller.Biophysics.BacterialModels':['CLBacterium.cl']},
)
