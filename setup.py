#!/usr/bin/env python

from distutils.core import setup

setup(name='CellModeller',
      packages=['CellModeller',
                'CellModeller.Biophysics',
                'CellModeller.Biophysics.BacterialModels',
                'CellModeller.Integration',
                'CellModeller.Regulation',
                'CellModeller.Signalling',
                'CellModeller.GUI'],
      package_data={'CellModeller.GUI':['test.ui', 'PyGLViewer.ui'],
                    'CellModeller.Integration':['CLCrankNicIntegrator.cl'],
                    'CellModeller.Integration':['CLEulerIntegrator.cl'],
                    'CellModeller.Biophysics.BacterialModels':['CLBacterium.cl']},
)
