#!/bin/bash

# plot pecube surface nodes history
# these files are in the output folder and start with "x_y_z_depth_"

for i in x_y_z_depth_*.dat; do
    infile=$i
    tempfile=${i/x_y_z_depth_0/Temps_plot_tec}
    pngfile=${infile/.dat/.png}

# only apply conversion if file does not exists
    if [[ ! -f $pngfile ]];
    then
        echo "infile: '$infile', pngfile: '$pngfile'"

        (gnuplot <<PLOT
        # save plot to png file, size: 1600 x 1200
        set terminal png size 1600,1200

        # set filename for png image file
        set output "$pngfile"

        # view from top, x and y coordinates
        set view map

        # set style for data file
        set style data points

        # set the size for all the data points
        set pointsize 0.9

        # show grid for x and y axis, line style -1 means black color
        set grid linestyle -1

        # set color model for color bar
        set palette rgbformulae 33,13,10

        # make two plots, two rows, one column
        set multiplot layout 2,1 title "files: $i, $tempfile

        set xlabel "x [km]"
        set xrange [0:950]
        set yrange [100:140]

        set title "moving surface nodes"
        set cbrange [1:900]
        set cblabel "id"
        plot "$infile" using 2:4:(\$3==0.0?\$1:1/0) palette

        set origin 0.0, 0.0
        set lmargin 0.0
        set rmargin 0.0

        set ylabel "z [km]"
        set xrange [-100:600]
        set yrange [-120:20]

        set title "temperature inside the model:"
        set cbrange [0:900]
        set cblabel "temperature [C]"
        splot "$tempfile" using 1:3:5 pointtype 7 palette notitle

        # set label for color boxes individually

#        set cblabel "x velocity [mm/year]"
#        set cbrange [-100:100]
#        splot "$tempfile" using 1:3:6 pointtype 7 palette notitle

#        set cblabel "z velocity [mm/year]"
#        set cbrange [-100:100]
#        splot "$tempfile" using 1:3:8 pointtype 7 palette notitle

PLOT
    ) &
    fi
done

