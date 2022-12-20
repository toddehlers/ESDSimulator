#!/bin/bash

# plot cascade tecplot files
# these files are in the output folder and start with "topo_tec_*.dat"

# 268.00000000    546.04026846 0.0       62150 0.0 0.0 0.0 0.0 0.0           0           0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0      1.54465004
#
# 0.24467583      0.96424759      0.00052581       62150 0.0 0.0 0.0 0.0 0.0           0           0 0.0 0.0 0.0      0.52581130 0.0 0.0      0.00052581 0.0 0.0 0.0 0.0 0.0      0.00039392

for i in topo_tec_*.dat; do
    infile=$i
    plotfile=${i/_tec/_plot_tec}
    pngfile=${plotfile/.dat/.png}

# only apply conversion if file does not exists
    if [[ ! -f $plotfile ]];
    then
        echo "infile: '$infile', plotfile: '$plotfile'"
# isolate interesting part of pecube output:
        awk 'NF == 24 {print $0}' $infile > $plotfile
    fi

# only call gnuplot if data file exists and png file does not exist
    if [[ ( -f $plotfile ) && ( ! -f $pngfile ) ]];
    then
        echo "data file: '$plotfile', png file: '$pngfile'"

        (gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        set view map

        # set interpolation
        set pm3d
        set dgrid3d 200,200,50

        # set labels for axis
        set xlabel "x [km]"
        set ylabel "y [km]"

        # set color model for color bar
        set palette rgbformulae 33,13,10

        # show grid for x and y axis, line style -1 means black color
        set grid linestyle -1

        set cblabel "z height [km]"

        set contour base

        set cntrparam levels 20

        # move legend + label outside
        set key outside

        splot "$plotfile" using 1:2:3 with lines nosurface
PLOT
    ) &
    fi


done

