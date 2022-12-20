#!/usr/bin/env python

import os
import sys
import fnmatch
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    print "usage: give base dir as command line parameter"
    sys.exit(1)

base_dir = sys.argv[1]

matches = []

for root, dirnames, filenames in os.walk(base_dir):
    for filename in fnmatch.filter(filenames, 'modtop25_50km.txt'):
        if "input" in root and "iter" in root:
            matches.append(os.path.join(root, filename))

matches.sort()

for file_name in matches:
    iteration = file_name.split("/")[1]
    z_values = []

    with open(file_name, "r") as input_file:
        z_values = [float(x) for x in input_file.readlines()]
        plt.plot(z_values, label=iteration)


plt.title("Changed topography (for each iteration)")
plt.xlabel("x [m]")
plt.ylabel("height z [m]")
plt.legend()
plt.show()

