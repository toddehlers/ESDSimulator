#!/bin/bash

counter=1

for folder in iter_0*; do
    echo "processing folder: $folder"
    error_plot="all_errors_${counter}.png"
    topo_plot="new_topo_${counter}.png"
    erosion_plot="all_erorion_${counter}.png"

    # echo $counter
    # echo $error_plot
    # echo $topo_plot
    # echo $erosion_plot

    gnuplot <<PLOT

    set key off
    set view map
    set terminal png size 1600,1200

    set output "$error_plot"
    splot "$folder/all_errors.txt" using 1:2:3 with points pointtype 5 pointsize 1 palette linewidth 30

    set output "$topo_plot"
    splot "$folder/new_topo.txt" every 100 using 1:2:3 with points pointtype 5 pointsize 0.5 palette linewidth 0.5

    set output "$erosion_plot"
    splot "$folder/all_erosion_rates.txt" using 1:2:3 with points pointtype 5 pointsize 1 palette linewidth 20
PLOT

    counter=$[$counter + 1]
done
