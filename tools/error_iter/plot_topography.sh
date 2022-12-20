#!/bin/bash

# Change these to your needs
TOPO_FILE="input/Pecube-D/oly_topo_v2.dat"

counter=1

for folder in iter_0*; do
    echo "processing folder: $folder"
    outputfile="topography_${counter}.png"

    awk 'BEGIN {i = 0; j = 0; k = 0;} {k = k + 1; if (k >= 100) {k = 0; print i, j, $1}; i = i + 1; if (i >= 7140) {i = 0; j = j + 1}}' $folder/$TOPO_FILE > $folder/topo.dat

    gnuplot <<PLOT
    set key off
    set view map
    set cbrange [-500.0:3000.0]
    set terminal png size 1600,1200

    set output "$outputfile"
    splot "$folder/topo.dat" using 1:2:3 with points pointtype 5 pointsize 0.5 palette linewidth 0.5
PLOT
    rm $folder/topo.dat
    counter=$[$counter + 1]
done