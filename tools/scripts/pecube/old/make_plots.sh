#!/bin/bash

#for i in data0*.dat; do
#    python plot_vec.py $i &
#done

#for i in grid_data0*.dat; do
#    python plot_vec.py $i &
#done

#for i in data0*.dat; do
#    gnuplot -e "inFile='$i'; outFile='$i.png'" velo_plot2.txt &
#done

for i in x_y_z_temp_*.dat; do
    gnuplot -e "inFile='$i'; outFile='$i.png'" temp_plot.txt &
done

# create animation afterwards:
# mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi

# rename files:
# for a in data_out_*.png; do b=$(printf data_out_%03d.png $(echo $a | tr -d [a-z_.])); mv $a $b; done

# find points in last file that are on the surface:
# awk '$5>60.0 {print $2;}' x_y_z_temp_0041.dat

# find points with specific point id:
# for i in x_y_z_*.dat; do awk '$2==8 {print $1, $2, $3, $5, $6}' $i; done > point_8.dat


