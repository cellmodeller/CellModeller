#!/bin/bash

THIS_DIR="`dirname $0`"

echo Installing CellModeller python package
$PYTHON setup.py install


# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
