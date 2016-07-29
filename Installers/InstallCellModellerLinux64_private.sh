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
echo 

# Change permission and run miniconda installer
echo
echo
echo 1. Installing Miniconda ...
echo ------------------------------------------------------------
echo
echo
chmod +x Miniconda-3.8.3-Linux-x86_64.sh
bash Miniconda-3.8.3-Linux-x86_64.sh -b

# Install dependencies
echo
echo
echo 2. Installing cellmodeller dependencies via miniconda
echo ------------------------------------------------------------
echo
echo
source ~/.bash_profile
source ~/.bashrc

#conda install pip
#pip install --upgrade pip
#pip install pyopengl
$HOME/miniconda/bin/conda install -yc https://conda.binstar.org/timrudge cellmodeller_dev

THIS_DIR="`dirname $0`"

echo
echo
echo 3. Setting up user environment
echo ------------------------------------------------------------
echo
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
source ~/.bashrc

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
