#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of IMEX-SfloW2d"

file="example6*"

rm -f $file

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dakota.txt"

if [ -f $file ] ; then
    rm $file
fi

