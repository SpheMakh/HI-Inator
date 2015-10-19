#!/bin/bash -ve

if [ -z "$USER" ]; then
  export USER=root
fi

echo "The configuratio file is: "
echo $CONFIG

pyxis CFG=${CONFIG} DESTDIR=/output OUTDIR=${OUTPUT} INDIR=${INPUT} azishe

