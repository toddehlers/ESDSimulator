#!/bin/bash

#for i in topo_tec_*.dat; do
#    tail -n +5 $i | head -n62150 > ${i/topo_tec_/xyz_} &
#done

for i in xyz_*.dat; do
    gnuplot -e "inFile='$i'; outFile='${i/dat/png}'" ../../plot1.txt &
done

