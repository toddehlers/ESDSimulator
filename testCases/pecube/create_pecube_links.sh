#!/bin/bash

echo "creating links for the pecube binary..."

for i in case_*;
do
    rm -v $i/pecube
    cp ../../pecube $i/pecube
    # ln -v -s /data/data2/willi/model_runs/ESD_Simulator/pecube $i/pecube
done

