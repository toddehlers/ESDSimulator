#!/bin/bash

# plot pecube velocity files
# these files are in the input folder and end with *.dat

for i in *_vel_new.dat; do
    infile=$i
    pngfile=${infile/.dat/.png}

# only apply conversion if file does not exists
    if [[ ! -f $pngfile ]];
    then
        echo "infile: '$infile', pngfile: '$pngfile'"

        gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        set view map

        # set labels for axis
        set xlabel "x [km]"
        set ylabel "z [km]"

        # set style for data file
        set style data points

        # set color model for color bar
        set palette rgbformulae 33,13,10

        # set new style for lines at index 100
        set style line 100 pointtype 7 palette

        # set the size for all the data points
        set pointsize 0.9

        # show grid for x and y axis, line style -1 means black color
        set grid linestyle -1

        # make two plots, two lines, one column
        set multiplot layout 2,1 title "velocity for $infile"

        # set label for color boxes individually

        set cblabel "x velocity [mm/year]"
        set cbrange [-60:60]
        splot "$infile" using 2:4:5 linestyle 100 notitle

        set cblabel "z velocity [mm/year]"
        set cbrange [-60:60]
        splot "$infile" using 2:4:7 linestyle 100 notitle
PLOT
    fi
done

