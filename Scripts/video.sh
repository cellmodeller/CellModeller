#!/bin/bash

# Note: this script requires imagemagick, ffmpeg and ghostscript
# The easiest way to install on OSX is to install homebrew package manager
# and then "brew install imagemagick" etc.

#for f in $( ls *.pickle ); do
for f in $( ls *.pdf ); do
    echo Processing: $f
#    cmpython ../../Scripts/Draw2DPDF.py $f
#    NAME=`basename $f .pickle`
    NAME=`basename $f .pdf`
#    convert $NAME-1um.pdf -resize 1024x1024 $NAME.png
    convert $NAME.pdf -resize 1024x1024 $NAME.png
done

ffmpeg -r 24 -i step-%04d0.png video.mp4
