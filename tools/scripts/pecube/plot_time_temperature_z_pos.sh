#!/bin/bash

awk '$1 == 1 {print $3, $4, $7}' time_temperature_history_0004.txt > plot_data.txt

input_file="plot_data.txt"

gnuplot <<PLOT

    # save plot to png file, size: 1600 x 1200
    set terminal png size 1600,1200

    # set filename for png image file
    set output "temperature_plot.png"

    # view from top, x and y coordinates
    set view map

    # set labels for axis
    set xlabel "time [Ma]"
    set ylabel "temperature [deg C]"

    # set style for data file
    set style data points

    set title "temperature plot"

    # activate second y axis
    set ytics nomirror
    set y2tics -25, 1, 5
    set y2label "depth [km]"

    plot "$input_file" using 1:2 title "temperature" axis x1y1, \
	 "$input_file" using 1:3 title "z position" axis x1y2

PLOT

