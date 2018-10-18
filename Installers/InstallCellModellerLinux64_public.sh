#!/bin/bash

clear
echo "********************************************************************************"
echo "*                                                                              *"
echo "*                          CellModeller Installer                              *"
echo "*                       ----------------------------                           *"
echo "*                                                                              *"
echo "*  This installer script will:                                                 *"
echo "*                                                                              *"
echo "*     1. Install CellModeller and dependencies into Anaconda python            *"
echo "*     2. Create your own ~/cellmodeller local CellModeller environment         *"
echo "*     3. Add some environment variables to your bash_profile/bashrc            *"
echo "*                                                                              *"
echo "********************************************************************************"
echo


# Install dependencies
echo
echo 1. Installing CellModeller and dependencies via Anaconda...
echo ------------------------------------------------------------
echo

command -v conda >/dev/null 2>&1 || { echo "You must install Anaconda before running this installer..." >&2; exit 1; }

conda create -yn cellmodeller 
source activate cellmodeller
conda install -yc trudge cellmodeller
source deactivate
source activate cellmodeller
conda install libgfortran==1

THIS_DIR="`dirname $0`"

echo
echo 2. Setting up user environment in $HOME/cellmodeller...
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
cp $THIS_DIR/PackageFiles/cmgui $BINDIR/
cp $THIS_DIR/PackageFiles/cmpython $BINDIR/

echo
echo 3. Adding paths to bash_profile
echo -------------------------------------------------------
echo
# Put an environment variable in bash_profile to tell CellModeller where things are
# Also append our bin dir to PATH
touch $HOME/.bashrc
touch $HOME/.bash_profile
echo export CMPATH=$CMDIR >> $HOME/.bash_profile
echo export "PATH=$BINDIR:$PATH" >> $HOME/.bash_profile
echo export CMPATH=$CMDIR >> $HOME/.bashrc
echo export "PATH=$BINDIR:$PATH" >> $HOME/.bashrc

