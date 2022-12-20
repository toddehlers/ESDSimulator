#!/bin/bash

# plot pecube velocity files
# these files are in the output folder and start with "Temps_tec"

for i in Temps_tec*.dat; do
    infile=$i
    outfile=${i/_tec/_plot_tec}
    pngfile=${outfile/.dat/.png}

# only apply conversion if file does not exists
    (
    if [[ ! -f $outfile ]];
    then
        echo "infile: '$infile', outfile: '$outfile'"
# isolate interesting part of pecube output:
        awk 'NF == 8 && $2 == 0.0 {print $0}' $infile > $outfile
    fi

# only call gnuplot if data file exists and png file does not exist
    if [[ ( -f $outfile ) && ( ! -f $pngfile ) ]];
    then
        echo "data file: '$outfile', png file: '$pngfile'"

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
        set multiplot layout 2,1 title "velocity for $outfile"

        # set label for color boxes individually

        set cblabel "x velocity [mm/year]"
        set cbrange [-100:100]
        splot "$outfile" using 1:3:6 linestyle 100 notitle

        set cblabel "z velocity [mm/year]"
        set cbrange [-100:100]
        splot "$outfile" using 1:3:8 linestyle 100 notitle
PLOT
    fi
    ) &
done

