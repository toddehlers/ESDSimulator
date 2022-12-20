#!/bin/bash

for i in 1 2 3 4 5 6 7 8 9; do
  echo $i
  python3 /esd/esd/data/willi/model_runs/ESDSimulator/tools/scripts/pecube/plot_time_temperature_history.py $i
  mv temperature_history_14.png temperature_history_14_p$i.png
done
