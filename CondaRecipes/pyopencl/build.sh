#!/bin/bash

git checkout tags/v2014.1
git submodule init
git submodule update
$PYTHON configure.py
$PYTHON setup.py build
make
$PYTHON setup.py install

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
