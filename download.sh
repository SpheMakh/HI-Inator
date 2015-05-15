#!/usr/bin/env bash

FILE=casapy-42.2.30986-1-64b.tar.gz
URL=https://svn.cv.nrao.edu/casa/linux_distro/old/$FILE

# make sure we are in the source folder
HERE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $HERE

if [ ! -e "$FILE" ]; then
    curl -O  $URL
else
    echo "$FILE already downloaded."
fi

# Some telescopes are not in the geodetic that comes with CASA
FILE=geodetic
URL=https://svn.cv.nrao.edu/svn/casa-data/distro/geodetic

# make sure we are in the source folder

if [ ! -d "$FILE" ]; then
    svn co  $URL
else
    echo "$FILE already downloaded."
fi
