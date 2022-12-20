#!/usr/bin/env python

# This python programm converts 2d move input files into files usable for Pecube, Matlab, Tecplot, gnuplot, etc.
# It takes user defined command line arguments and exports the resulting topography
#
# Written by Willi Kappler, willi.kappler@uni-tuebingen.de
#

import sys
import getopt
import time
import math
import inspect

VERSION = "V0.4 (2014.03.14)"

def usage(errorCode):
    executable = sys.argv[0]
    print
    print "makeTopo {0:s}".format(VERSION)
    print "usage: {0:s}".format(executable)
    print "-d deformation front in km (x axis)"
    print "-e elevation left of def. front in km (z axis)"
    print "-a topo angle in degree"
    print "-i ID for line in 2d move file (integer value / list of integer values)"
    print "-s step size in x direction for output file (in km)"
    print "-x maximum (last) x value on the right side (in km)"
    print "-f filename of 2d move file"
    print "-o output filename"
    print "-h print this information"
    print "-b x position of first point"
    print
    print "examples:"
    print "{0:s} -d 100 -e 5 -a 2 -i 12 -f topo1.dat -o surface1.dat".format(executable)
    print
    print "{0:s} -d 10 -e 0 -a 2.3 -i 51,55,73,79 -f topo2.dat -o surface2.dat".format(executable)
    print
    print "the 2D Move input file has the following layout:"
    print "a comment starts with '#'"
    print "4 columns are separated by spaces in that order: x, y, z, ID"

    sys.exit(errorCode)

# This functions processes the input data (list of tuples) and returns the result as a list of tuples
def processData(inputValues, deformationFront, elevation, topoAngle, inputFilename, outputFilename, xStep, xMax, startX):
    result = []

    # initialize current point
    px = startX
    py = 0.0
    pz = 0.0

    # sanity check
    if deformationFront < px:
        print "deformation front smaller than the first point: {0:f} < {1:f}".format(deformationFront, px)
        sys,exit(1)

    # indices of values on the left and on the right of the current point
    leftIndex = 0
    rightIndex = 1

    zFactor = math.tan(topoAngle * math.pi / 180.0)

    crossPointReached = False
    lastDiff = 0.0
    currentDiff = 0.0

    # main loops:
    # here we go from left to right (along the x axis) and calculate the correct z value for each point

    # first the deformation front
    while px <= deformationFront:
        result.append((px, py, elevation))
        px = px + xStep

    # x values of input: inputValues[_][0]
    # z values of input: inputValues[_][2]

    # find index position in input data, so that current px matches it
    while px > inputValues[rightIndex][0]:
        rightIndex = rightIndex + 1

    leftIndex = rightIndex - 1

    print "after deformation front:"
    print "rightIndex: {0:d}, leftIndex: {1:d}, {2:f}, {3:f}".format(rightIndex, leftIndex, inputValues[rightIndex][0], px)

    # now everything that is to the left of the deformation front
    while px <= inputValues[-1][0]:
        # calculate linear interpolated z value
        if (inputValues[rightIndex][0] - inputValues[leftIndex][0]) == 0:
          print "input values invalid, the same x position appeas several times: {0:f}".format(inputValues[rightIndex][0])
          print "px: {0:f}".format(px)
          sys.exit(1)
        else:
          pzInter = inputValues[leftIndex][2] + ((inputValues[rightIndex][2] - inputValues[leftIndex][2]) * (px - inputValues[leftIndex][0]) / (inputValues[rightIndex][0] - inputValues[leftIndex][0]))

        # calculate z value with the given angle
        pz = elevation + ((px - deformationFront) * zFactor)

        if pz > pzInter:
            if pzInter > 0.0:
                pz = pzInter
            else:
                pz = 0.0

        result.append((px, py, pz))

        px = px + xStep

        # move to next pair of points in input file
        if px > inputValues[rightIndex][0]:
            leftIndex = leftIndex + 1
            rightIndex = rightIndex + 1

    # here we fill the right side to fit into the pecube model
    while px <= xMax:
        # we take the last z value
        result.append((px, py, pz))

        px = px + xStep


    return result

def main(argv):
    deformationFront = None
    elevation = None
    topoAngle = None
    lineIDs = None
    xStep = None
    xMax = None
    inputFilename = None
    outputFilename = None
    startX = 0.0

    print "{0:s}\n{1:s}".format(__file__, inspect.getfile(inspect.currentframe()))
    print "version: {0:s}".format(VERSION)

    # parse all command line arguments here
    try:
        opts, args = getopt.getopt(argv, "d:e:a:i:f:s:x:o:h:b:", [])
        for opt, arg in opts:
            if opt == "-d":
                deformationFront = float(arg)
                print "deformation front: {0:f}".format(deformationFront)
            elif opt == "-e":
                elevation = float(arg)
                print "elevation: {0:f}".format(elevation)
            elif opt == "-a":
                topoAngle = float(arg)
                print "topoAngle: {0:f}".format(topoAngle)
            elif opt == "-i":
                # check if line ID is a single integer value or a list of integer values
                if "," in arg:
                    lineIDs = [int(newID) for newID in arg.split(",")]
                else:
                    lineIDs = [int(arg)]
                print "line id: {0:s}".format(str(lineIDs))
            elif opt == "-f":
                inputFilename = arg
                print "input filename: {0:s}".format(inputFilename)
            elif opt == "-s":
                xStep = float(arg)
                print "x step: {0:f}".format(xStep)
            elif opt == "-x":
                xMax = float(arg)
                print "max x: {0:f}".format(xMax)
            elif opt == "-o":
                outputFilename = arg
                print "output filename: {0:s}".format(outputFilename)
            elif opt == "-h":
                usage(0)
            elif opt == "-b":
                startX = float(arg)
                print "x position of first point"
    except getopt.GetoptError as e:
        print "unknown command line argument:"
        print e
        usage(1)
    except ValueError as e:
        print "invalid command line argument:"
        print e
        usage(1)

    # check if all command line arguments are given
    if deformationFront == None:
        print "deformation front missing!"
        usage(1)
    if elevation == None:
        print "elevation missing"
        usage(1)
    if topoAngle == None:
        print "topo angle missing"
        usage(1)
    if lineIDs == None:
        print "line ID missing"
        usage(1)
    if xStep == None:
        print "x step missing"
        usage(1)
    if xMax == None:
        print "max x missing"
        usage(1)
    if inputFilename == None:
        print "2d move input filename missing"
        usage(1)
    if outputFilename == None:
        print "output filename missing"
        usage(1)

    # open input file and read in data

    lineValues = {}
    unchangedValues = []

    for lID in lineIDs:
        lineValues[lID] = []

    try:
        inputFile = open(inputFilename)

        for line in inputFile.readlines():
            line = line.strip()
            # ignore comments and blank lines
            if line.startswith("#") or len(line) == 0:
                continue
            else:
                try:
                    (px, py, pz, ID) = line.split()
                    newID = int(ID)
                    if newID in lineIDs:
                        lineValues[newID].append((float(px), float(py), float(pz)))
                    else:
                        unchangedValues.append(line)
                except ValueError as e:
                    print e
                    print
                    print "syntax error in line:"
                    print "\n"
                    print line
                    print "\n"
                    print "the format should be px, py, pz, ID"
                    print "use '#' to comment out lines line this:"
                    print "# This is a comment. Exported data from 2D Move"
        inputFile.close()
    except IOError as e:
        print "could not read input file '{0:s}'".format(inputFilename)
        print e
        sys.exit(1)

    for lID in lineIDs:
        # check if we have enough points / if the line id was correct
        if len(lineValues[lID]) < 3:
            print "too few points in input file for line ID: '{0:d}', maybe line ID is not correct ?".format(lID)
            sys.exit(1)

        # sort line by inc. x values
        lineValues[lID].sort()

    # process the line values
    result = {}

    for lID in lineIDs:
        result[lID] = processData(lineValues[lID], deformationFront, elevation, topoAngle, inputFilename, outputFilename, xStep, xMax, startX)

    # write result to output file:
    try:
        outputFile = open(outputFilename, "w")

        for lID in lineIDs:
            for (px, py, pz) in result[lID]:
                outputFile.write("{0:f}, {1:f}, {2:f} {3:d}\n".format(px, py, pz, lID))

        # only write unchanged lines if the user provided more than one line id
        if len(lineIDs) > 1:
            for line in unchangedValues:
                outputFile.write(line + "\n")

        outputFile.close()
    except IOError as e:
        print "could not open output file '{0:s}'".format(outputFilename)
        print e
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])

