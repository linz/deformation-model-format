#!/bin/sh
set -e
srcdir='.'
tooldir='../../../tools'
if [ -n "$1" ]; then
    srcdir=$1
fi
if [ -n "$2" ]; then
    tooldir=$2
fi
for f in $srcdir/*.tif.json; do
    t=`basename $f .json`
    echo "Building $t"
    $tooldir/deformation_csv_to_gtiff.py $f $srcdir/$t
done
