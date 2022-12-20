#!/bin/bash

# calculate pecube input surface files min and max values
# these files are in the input folder and usually start with "topo_"

for i in topo_*.dat; do
    echo "values for $i:"
    awk 'BEGIN {m1 = -100000.0; m2 = 100000.0}; {if ($1 > m1) m1 = $1; if ($1 < m2) m2 = $1}; END {printf("xmin: %f, xmax: %f\n", m2, m1);}' $i
    awk 'BEGIN {m1 = -100000.0; m2 = 100000.0}; {if ($3 > m1) m1 = $3; if ($3 < m2) m2 = $3}; END {printf("zmin: %f, zmax: %f\n", m2, m1);}' $i
done
