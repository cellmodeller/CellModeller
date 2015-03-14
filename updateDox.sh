#!/bin/bash

# Run doxygen
echo 1. Run doxygen on cellmodeller source code...
doxygen cellmodellerDox

# Bring html output to remote gh-pages repo
echo --
echo 2. Bring new doxygen html output to remote gh-page
git add html
git add updateDox.sh
git commit -m "update html on master"
git checkout -b gh-pages origin/gh-pages
git checkout master -- html
git add html
git commit -m "update /html on gh-pages"
git push
git checkout master
git branch -d gh-pages


