#!/bin/bash

THIS_DIR="`dirname $0`"

# How to get release version from git tag?
VERSION=`git describe`

# OSX from public repo
INSTALLERDIR=CellModeller-$VERSION-OSX
mkdir -p $INSTALLERDIR
mkdir -p $INSTALLERDIR/PackageFiles
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Scripts $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Doc $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmpython $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmgui $INSTALLERDIR/PackageFiles/
cp -f $THIS_DIR/InstallCellModellerOSX_public.command $INSTALLERDIR/InstallCellModeller
cp -f $THIS_DIR/UninstallCellModeller $INSTALLERDIR/UninstallCellModeller
chmod +x $INSTALLERDIR/InstallCellModeller
tar -czf $INSTALLERDIR.tgz $INSTALLERDIR 
rm -rf $INSTALLERDIR

# OSX from private repo
INSTALLERDIR=CellModeller-$VERSION-OSX-dev
mkdir -p $INSTALLERDIR
mkdir -p $INSTALLERDIR/PackageFiles
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Scripts $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Doc $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmpython $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmgui $INSTALLERDIR/PackageFiles/
cp -f $THIS_DIR/InstallCellModellerOSX_private.command $INSTALLERDIR/InstallCellModeller
cp -f $THIS_DIR/UninstallCellModeller $INSTALLERDIR/UninstallCellModeller
chmod +x $INSTALLERDIR/InstallCellModeller
tar -czf $INSTALLERDIR.tgz $INSTALLERDIR 
rm -rf $INSTALLERDIR

# Linux 64bit from public repo
INSTALLERDIR=CellModeller-$VERSION-Linux64
mkdir -p $INSTALLERDIR
mkdir -p $INSTALLERDIR/PackageFiles
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Scripts $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Doc $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmpython $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmgui $INSTALLERDIR/PackageFiles/
cp -f $THIS_DIR/InstallCellModellerLinux64_public.sh $INSTALLERDIR/InstallCellModeller.sh
cp -f $THIS_DIR/UninstallCellModeller $INSTALLERDIR/UninstallCellModeller
chmod +x $INSTALLERDIR/InstallCellModeller.sh
tar -czf $INSTALLERDIR.tgz $INSTALLERDIR 
rm -rf $INSTALLERDIR

# Linux 64bit from private repo
INSTALLERDIR=CellModeller-$VERSION-Linux64-dev
mkdir -p $INSTALLERDIR
mkdir -p $INSTALLERDIR/PackageFiles
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Scripts $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Examples $INSTALLERDIR/PackageFiles/
cp -Rf $THIS_DIR/../Doc $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmpython $INSTALLERDIR/PackageFiles/
cp $THIS_DIR/cmgui $INSTALLERDIR/PackageFiles/
cp -f $THIS_DIR/InstallCellModellerLinux64_private.sh $INSTALLERDIR/InstallCellModeller.sh
cp -f $THIS_DIR/UninstallCellModeller $INSTALLERDIR/UninstallCellModeller
chmod +x $INSTALLERDIR/InstallCellModeller.sh
tar -czf $INSTALLERDIR.tgz $INSTALLERDIR 
rm -rf $INSTALLERDIR

