#!/bin/bash

# plot pecube ages files
# these files are in the output folder and start with "Ages_tec"

for i in Ages_tec*.dat; do
    infile=$i
    outfile=${i/_tec/_plot_tec}
    pngfile=${outfile/.dat/.png}

# only apply conversion if file does not exists
    if [[ ! -f $outfile ]];
    then
        echo "infile: '$infile', outfile: '$outfile'"
# isolate interesting part of pecube output:
        awk 'NF > 4 && $2 == 0.0 {print $0}' $infile > $outfile
    fi

# only call gnuplot if data file exists and png file does not exist
    if [[ ( -f $outfile ) && ( ! -f $pngfile ) ]];
    then
        echo "data file: '$outfile', png file: '$pngfile'"

        (gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        set view map

        # set labels for axis
        set xlabel "x [km]"
        set ylabel "age [Ma]"

        # set style for data file
        set style data points

        set yrange [-20:150]

        set title "file: $infile"

        plot "$outfile" using 1:4 title "AHe Age (Ma)", \
             "$outfile" using 1:5 title "AFT Age (Ma)", \
             "$outfile" using 1:6 title "ZHe Age (Ma)", \
             "$outfile" using 1:7 title "ZFT Age (Ma)", \
             "$outfile" using 1:8 title "MAr Age (Ma)"
PLOT
    ) &
    fi
done

