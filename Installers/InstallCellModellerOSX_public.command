#!/bin/bash

clear
echo "********************************************************************************"
echo "*                                                                              *"
echo "*                          CellModeller Installer                              *"
echo "*                       ----------------------------                           *"
echo "*                                                                              *"
echo "*  This installer script will:                                                 *"
echo "*                                                                              *"
echo "*     1. Install miniconda into ~/miniconda                                    *"
echo "*     2. Install CellModeller and dependencies into this miniconda python      *"
echo "*     3. Create your own ~/cellmodeller local CellModeller environment         *"
echo "*     4. Add some environment variables to your bash_profile/bashrc            *"
echo "*                                                                              *"
echo "********************************************************************************"
echo

# Download miniconda
echo
echo 1. Downloading Miniconda...
echo ---------------------------
echo
curl https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh > Miniconda.sh

# Change permission and run miniconda installer
echo
echo 2. Installing Miniconda ...
echo ---------------------------
echo
chmod +x Miniconda.sh
bash Miniconda.sh -b
rm Miniconda.sh

# Install dependencies
echo
echo 3. Installing CellModeller and dependencies via miniconda...
echo ------------------------------------------------------------
echo
source ~/.bash_profile

$HOME/miniconda/bin/conda install -yc https://conda.binstar.org/timrudge cellmodeller

THIS_DIR="`dirname $0`"

echo
echo 4. Setting up user environment in $HOME/cellmodeller...
echo -------------------------------------------------------
echo
# Setup the users cellmodeller directory with scripts etc.
#
# TODO: user option here to choose directory location
#
CMDIR=$HOME/cellmodeller
DATADIR=$CMDIR/data
MODELDIR=$CMDIR/Models
EXAMPLEDIR=$CMDIR/Examples
SCRIPTDIR=$CMDIR/Scripts
DOCDIR=$CMDIR/Doc
BINDIR=$CMDIR/bin

# Create CellModeller dir in user home area
mkdir -p $CMDIR
mkdir -p $DATADIR
mkdir -p $BINDIR
mkdir -p $MODELDIR
#mkdir -p $EXAMPLEDIR
#mkdir -p $SCRIPTDIR
#mkdir -p $DOCDIR

# Copy required stuff
cp -Ri $THIS_DIR/PackageFiles/Scripts $CMDIR
cp -Ri $THIS_DIR/PackageFiles/Doc $CMDIR
cp -Ri $THIS_DIR/PackageFiles/Examples $CMDIR
cp $THIS_DIR/PackageFiles/cmpython $BINDIR/
cp $THIS_DIR/PackageFiles/cmgui $BINDIR/

# Put an environment variable in bash_profile to tell CellModeller where things are
# Also append our bin dir to PATH
echo export CMPATH=$CMDIR >> $HOME/.bash_profile
echo export PATH=$BINDIR:$PATH >> $HOME/.bash_profile
source ~/.bash_profile

clear
echo 
echo 
echo "********************************************************************************"
echo "*                                                                              *"
echo "*                           SUCCESS! All done.                                 *"
echo "*                                                                              *"
echo "********************************************************************************"
echo 
echo 
