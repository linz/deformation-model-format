#!/bin/sh
for f in *.tif.json; do
    t=`basename $f .json`
    echo "Building $t"
    ../../../tools/deformation_csv_to_gtiff.py $f ../$t
done
