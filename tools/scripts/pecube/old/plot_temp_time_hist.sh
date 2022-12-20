#!/bin/bash

# plot pecube temperature over time history for the nodes on the surface
# these files are in the output folder and start with "Plot_hist_"
# you need to create them with the python script "plot_temp_time_history.py"

for i in Plot_hist_*.dat; do
    infile=$i
    pngfile=${infile/.dat/.png}

# only call gnuplot if data file exists and png file does not exist
    if [[ ! -f $pngfile ]];
    then
        echo "data file: '$infile', png file: '$pngfile'"

        (gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        # set view map

        # set labels for axis
        set xlabel "time [Ma]"
        set ylabel "temperature [C]"

        # set style for data file
        set style data points

        set xrange [-50:0]

        set xtics 2

        set mxtics
        set mytics

        # make two plots, two lines, one column
        set multiplot layout 2,1 title "infile $infile"

        set ylabel "temperature [C]"
        set yrange [0:400]
        set ytics 20
        plot "$infile" using 1:2, "$infile" using 1:3, "$infile" using 1:4, "$infile" using 1:5

        set ylabel "z [km]"
        set yrange [80:140]
        set ytics 5
        plot "$infile" using 1:6, "$infile" using 1:7, "$infile" using 1:8, "$infile" using 1:9
PLOT
    ) &
    fi
done
