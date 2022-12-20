#!/bin/bash

# plot pecube surface files
# these files are in the output folder and start with "surface_temp_"
# also plot the temperaturs from the file Temps_tec*

for i in surface_temp_*.dat; do
    pngfile=${i/.dat/.png}
    datafile2=${i/surface_temp_0/Temps_plot_tec}
    agesfile=${i/surface_temp_0/Ages_plot_tec}
    currentStep=${i:14:3}

    if [[ ! -f $pngfile ]];
    then
        echo "data file: '$i', data file2: '$datafile2', ages file: '$agesfile', png file: '$pngfile'"

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
        set xrange [-100:600]

        set key outside right

        # make four plots, two rows, one column
        set multiplot layout 2,1 title "files: $i, $datafile2, $agesfile"

        set title "ages:"
        set ylabel "age [Ma]"
        set yrange [-20:150]
        set lmargin 10.0
        set rmargin 21.5
        plot "$agesfile" using 1:4 title "AHe Age (Ma)", \\
             "$agesfile" using 1:5 title "AFT Age (Ma)", \\
             "$agesfile" using 1:6 title "ZHe Age (Ma)", \\
             "$agesfile" using 1:7 title "ZFT Age (Ma)", \\
             "$agesfile" using 1:8 title "MAr Age (Ma)"

#        set title "temperature on the surface:"
#        set ylabel "temperature [C]"
#        set yrange [0:500]
#        plot "$i" using (\$1 == $currentStep + 1 ? (\$3 - 65) : 1/0):(\$4 == 0.0 ? \$5 : 1/0) title "step $(echo $currentStep + 1 | bc)", \\
#             "$i" using (\$1 == $currentStep + 2 ? (\$3 - 65) : 1/0):(\$4 == 0.0 ? \$5 : 1/0) title "step $(echo $currentStep + 2 | bc)", \\
#             "$i" using (\$1 == $currentStep + 3 ? (\$3 - 65) : 1/0):(\$4 == 0.0 ? \$5 : 1/0) title "step $(echo $currentStep + 3 | bc)"
#
#        set title "temperature at different depth:"
#        set ylabel "temperature [C]"
#        plot "$datafile2" using 1:(\$3 == -1.026 ? \$5 : 1/0) title "-1.026 km", \\
#             "$datafile2" using 1:(\$3 == -2.051 ? \$5 : 1/0) title "-2.051 km", \\
#             "$datafile2" using 1:(\$3 == -3.077 ? \$5 : 1/0) title "-3.077 km", \\
#             "$datafile2" using 1:(\$3 == -4.103 ? \$5 : 1/0) title "-4.103 km", \\
#             "$datafile2" using 1:(\$3 == -5.128 ? \$5 : 1/0) title "-5.128 km"

        set title "temperature inside the model:"
        set ylabel "z [km]"
        set yrange [-120:20]
        set cbrange [0:900]
        set cblabel "temperature [C]"
        set origin 0.0, 0.0
        set lmargin 0.0
        set rmargin 0.0
        splot "$datafile2" using 1:3:5 pointtype 7 palette notitle

PLOT
    ) &
    fi
done

