#!/usr/bin/env python

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

def processFile(fileName):
    xVals = []
    yVals = []
    zVals = []

    with open(fileName, "r") as inputFile:
        for line in inputFile:
            lineData = line.split()
            if len(lineData) == 24:
                x = float(lineData[0])
                y = float(lineData[1])
                z = float(lineData[2])

                if (z < -5.0):
                    z = 5.0
                if (z > 5.0):
                    z = 5.0

                xVals.append(x)
                yVals.append(y)
                zVals.append(z)

    return (xVals, yVals, zVals)

def plotData(fileName, xVals, yVals, zVals):
     plt.figure(1).clear()
     plt.figure(1).subplots_adjust(top=0.96, bottom=0.06, left=0.06, right=0.98)

     ax = plt.subplot(111)

     numOfContours = 15
     contourFactor = 6.5 / float(numOfContours - 3)

     plt.tricontourf(xVals, yVals, zVals, [((float(x) * contourFactor) - 5.0) for x in range(numOfContours)])

     cb = plt.colorbar(ticks=[-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0])

     plt.tricontour(xVals, yVals, zVals, [((float(x) * contourFactor) - 5.0) for x in range(numOfContours)], colors="k")

     plt.xlabel("x pos [km]")
     plt.ylabel("y pos [km]")
     plt.axis([-100, 2300, -100, 2100])

     ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2)

     ax.xaxis.set_major_locator(mti.MultipleLocator(200))
     ax.xaxis.set_minor_locator(mti.MultipleLocator(50))

     ax.yaxis.set_major_locator(mti.MultipleLocator(100))
     ax.yaxis.set_minor_locator(mti.MultipleLocator(50))

     ax.set_title("input file: {0:s}".format(fileName))

     baseFileName = os.path.basename(fileName)
     fNameWithoutExt = os.path.splitext(baseFileName)[0]

     print "saving image file: '{0:s}.png'".format(fNameWithoutExt)

     # plt.show()

     plt.savefig("{0:s}.png".format(fNameWithoutExt), dpi=300)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print "usage: "
        print "{0}s input_file_name".format(sys.argv[0])
        sys.exit(1)

    # global settings for matplotlib
    mpl.rcParams["font.size"] = 6
    mpl.rcParams["lines.markersize"] = 3.0
    mpl.rcParams["lines.linewidth"] = 0.25
    mpl.rcParams["contour.negative_linestyle"] = "solid"

    for fileName in sys.argv[1:]:
        print "processing file: '{0:s}'".format(fileName)

        xVals, yVals, zVals = processFile(fileName)

        print "num of data values: x: {0:d}, y: {1:d}, z: {2:d}".format(len(xVals), len(yVals), len(zVals))
        print "min(x): {0:f}, max(x): {1:f}".format(min(xVals), max(xVals))
        print "min(y): {0:f}, max(y): {1:f}".format(min(yVals), max(yVals))
        print "min(z): {0:f}, max(z): {1:f}".format(min(zVals), max(zVals))

        plotData(fileName, xVals, yVals, zVals)
