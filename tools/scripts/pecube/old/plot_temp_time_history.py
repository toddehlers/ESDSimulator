#!/usr/bin/env python

import glob
import os.path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

timeFileNames = glob.glob("./Time_*.dat")

# columns = [10, 101, 1621, 2967]
columns = [10, 51, 74, 103, 127, 153, 182, 206]
# columns = [60, 65, 70, 71, 72, 73, 74]

for timeFileName in timeFileNames:
    tempFileName = timeFileName.replace("Time_", "Temp_", 1)
    zFileName = timeFileName.replace("Time_", "Z_", 1)
    xFileName = timeFileName.replace("Time_", "X_", 1)

    print "checking files:\n {0:s},\n {1:s},\n {2:s},\n {3:s}".format(timeFileName, tempFileName, zFileName, xFileName)

    if not os.path.isfile(timeFileName):
        print "could not open file: {0:s}".format(timeFileName)
        sys.exit(1)

    if not os.path.isfile(tempFileName):
        print "could not open file: {0:s}".format(tempFileName)
        sys.exit(1)

    if not os.path.isfile(zFileName):
        print "could not open file: {0:s}".format(zFileName)
        sys.exit(1)

    if not os.path.isfile(xFileName):
        print "could not open file: {0:s}".format(xFileName)
        sys.exit(1)

    timeValues = []
    tempValues = []
    zValues = []
    xValues = []

    with open(timeFileName, "r") as timeFile:
        with open(tempFileName, "r") as tempFile:
            with open(zFileName, "r") as zFile:
                with open(xFileName, "r") as xFile:

                    for lineTimeFile in timeFile:
                        timeValues.append(- float(lineTimeFile.split()[0]))

                    for lineTempFile in tempFile:
                        line = lineTempFile.split()
                        values = []
                        for i in columns:
                            values.append(line[i])
                        tempValues.append(values)

                    for lineZFile in zFile:
                        line = lineZFile.split()
                        values = []
                        for i in columns:
                            values.append(line[i])
                        zValues.append(values)

                    for lineXFile in xFile:
                        line = lineXFile.split()
                        values = []
                        for i in columns:
                            values.append(line[i])
                        xValues.append(values)

    if len(timeValues) != len(tempValues):
        print "number of values are not equal for files: {0:s} and {1:s}".format(timeFileName, tempFileName)
        sys.exit(1)

    if len(timeValues) != len(zValues):
        print "number of values are not equal for files: {0:s} and {1:s}".format(timeFileName, zFileName)
        sys.exit(1)

    if len(timeValues) != len(xValues):
        print "number of values are not equal for files: {0:s} and {1:s}".format(timeFileName, xFileName)
        sys.exit(1)

    # global setting for all plots
    # see http://matplotlib.org/users/customizing.html
    mpl.rcParams["font.size"] = 6
    mpl.rcParams["figure.figsize"] = 5, 8 # width and height
    mpl.rcParams["lines.markersize"] = 3.0
    mpl.rcParams['axes.color_cycle'] = ["#0000ff", "#00aa00", "#aa0000", "#00aaaa", "#aa00aa", "#aaaa00", "#000000", "#8888ff"]

    plt.figure(1).clear()
    plt.figure(1).subplots_adjust(top=0.96, bottom=0.04, left=0.08, right=0.98)

    ax = plt.subplot(311)

    for i in xrange(len(columns)):
        plt.plot(timeValues, [items[i] for items in tempValues], "+", label="column: {0:d}".format(columns[i]))

    plt.xlabel("time [Ma]")
    plt.ylabel("temperature [C]")
    plt.axis([-50, 0, 0, 400])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4)
    ax.xaxis.set_minor_locator(mti.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mti.MultipleLocator(10))
    ax.set_title("input file: {0:s}".format(timeFileName))

    ax = plt.subplot(312)

    for i in xrange(len(columns)):
        plt.plot(timeValues, [items[i] for items in zValues], "+", label="column: {0:d}".format(columns[i]))

    plt.xlabel("time [Ma]")
    plt.ylabel("z [km]")
    plt.axis([-50, 0, 80, 140])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4)
    ax.xaxis.set_minor_locator(mti.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mti.MultipleLocator(2))

    ax = plt.subplot(313)

    for i in xrange(len(columns)):
        plt.plot([items[i] for items in xValues], [items[i] for items in zValues], "+", label="column: {0:d}".format(columns[i]))

    plt.xlabel("x [km]")
    plt.ylabel("z [km]")
    plt.axis([0, 700, 80, 140])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4)
    ax.xaxis.set_minor_locator(mti.MultipleLocator(10))
    ax.yaxis.set_minor_locator(mti.MultipleLocator(2))

    outFileName = timeFileName.replace("Time_", "Plot_", 1)
    plt.savefig(outFileName.replace(".dat", ".png", 1), dpi=200)

#    break
