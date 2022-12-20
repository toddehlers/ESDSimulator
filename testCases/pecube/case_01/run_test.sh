#!/bin/bash

CURRENT_DATE=$(date +"%Y_%m_%d__%H_%M")

./pecube %> out_$CURRENT_DATE.txt

mkdir output/Pecube_$CURRENT_DATE
mv output/Pecube\-D/* output/Pecube_$CURRENT_DATE

