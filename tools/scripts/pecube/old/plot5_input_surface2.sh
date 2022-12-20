#!/bin/bash

# plot pecube input surface files
# these files are in the input folder and usually start with "topo_"

for i in interpol_surface_*.dat; do
    pngfile=${i/.dat/.png}

    if [[ ! -f $pngfile ]];
    then
        echo "data file: '$i', png file: '$pngfile'"

        (gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        set view map

        # set labels for axis
        set xlabel "x [km]"

        # set style for data file
        set style data points

        # set the size for all the data points
        set pointsize 0.9

        # show grid for x and y axis, line style -1 means black color
        set grid linestyle -1

        # set color model for color bar
        set palette rgbformulae 33,13,10

        # set x range for both plots
        set xrange [-100:800]

        set key outside bottom center

        set title "surface"
        set ylabel "z [km]"
        set yrange [-10:10]

        plot "$i" using 1:3

PLOT
    ) &
    fi
done

