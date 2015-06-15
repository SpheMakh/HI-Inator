#!/bin/bash -ve

if [ -z "$1" ]; then
    DATA=/
else
    DATA=$1
fi

if [ -z "$USER" ]; then
  export USER=root
fi

echo "The configuratio file is: "
echo $config

pyxis CFG=${DATA}/input/${config} DESTDIR=${DATA}/output OUTDIR=${DATA}/output azishe

