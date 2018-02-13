#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of IMEX-SfloW2d"

file="example4*"

rm -f $file

file="IMEX_SfloW2D.inp"

if [ -f $file ] ; then
    rm $file
fi

file="topography_dem.asc"

if [ -f $file ] ; then
    rm $file
fi

file="pile.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

