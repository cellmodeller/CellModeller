#!/bin/bash

# Note: this script requires imagemagick, ffmpeg and ghostscript
# The easiest way to install on OSX is to install homebrew package manager
# and then "brew install imagemagick" etc.

# Run Draw2DPDF to generate pdf files
for f in $( ls *.pickle ); do
    echo Processing: $f
    python $HOME/Code/CellModeller/Scripts/Draw2DPDF.py $f
done

# Convert and resize etc. pdf files into jpegs
for f in $( ls *.pdf ); do
    NAME=`basename $f .pdf`
    convert \
           -colorspace RGB \
	   -verbose       \
           -density 150   \
            $NAME.pdf      \
            $NAME.png
done

# Run ffmpeg to generate video file
ffmpeg -framerate 7 -i %*.png -vf scale=1920:1080 -r 24 $1
