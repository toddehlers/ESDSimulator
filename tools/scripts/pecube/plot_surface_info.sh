#!/bin/bash

for ID in 2012 2013 2014 2015 2016 2017 2018 2062 2063; do
  echo "processind id: $ID"
  awk -v id="$ID" '$2 == id' output/Pecube-D/surface_information.dat > id_$ID.txt
  python3 /esd/esd/data/willi/model_runs/ESDSimulator/tools/scripts/pecube/plot_surface_info.py id_$ID.txt
done
