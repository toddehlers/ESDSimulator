#!/usr/bin/env python

import os.path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

# read in tecplot file
xpos = []
AHe = []
ZHe = []
AFT = []
ZFT = []

if len(sys.argv) != 2:
    print "usage:"
    print "{0:s input_filenmae}".format(sys.argv[0])
    sys.exit(1)

tecplotFileName = sys.argv[1]

with open(tecplotFileName, "r") as tecplotFile:
    lineCounter = 0
    for line in tecplotFile:
        lineCounter = lineCounter + 1
        if lineCounter > 3:
            values = line.split()
            xpos.append(values[0])
            AHe.append(values[3])
            ZHe.append(values[4])
            AFT.append(values[5])
            ZFT.append(values[6])

# global setting for all plots
# see http://matplotlib.org/users/customizing.html
mpl.rcParams["font.size"] = 6
mpl.rcParams["lines.markersize"] = 3.0
mpl.rcParams['axes.color_cycle'] = ["#0000ff", "#00aa00", "#aa0000", "#00aaaa", "#aa00aa", "#aaaa00", "#000000", "#8888ff"]

plt.figure(1).clear()
plt.figure(1).subplots_adjust(top=0.96, bottom=0.06, left=0.06, right=0.98)

ax = plt.subplot(111)

plt.plot(xpos, AHe, "+", label="AHe")
plt.plot(xpos, ZHe, "+", label="ZHe")
plt.plot(xpos, AFT, "+", label="AFT")
plt.plot(xpos, ZFT, "+", label="ZFT")

plt.xlabel("x pos [km]")
plt.ylabel("age [Ma]")
plt.axis([0, 650, 60, 110])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2)
ax.xaxis.set_minor_locator(mti.MultipleLocator(5))
ax.yaxis.set_minor_locator(mti.MultipleLocator(1))
ax.set_title("input file: {0:s}".format(tecplotFileName))

plt.savefig("ages.png", dpi=200)

