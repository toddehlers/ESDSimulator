#!/usr/bin/env python

import sys
import os

# velocity.py
# Written by Willi Kappler
#
# This program reads an input file with the following format:
#
# first line: depth in y dimension in km
# second line: number of nodes in y dimension
# third line: number of samples to take for interpolation
# fourth line: step size in x direction for interpolation
# fifth line to end of file: time in Myrs and filename
#
# example:
#
# 10.0
# 5
# 10
# 1.0
# 50.0 file1.dat
# 25.0 file2.dat
# 12.56 file3.dat
# 7.1 file4.dat
# 1.8 file5.dat
# 0.0 file6.dat
#
# The input files specified are output files from 2d Move
# Those files must have the following format for each line:
#
# x y z color id
#
# The output file will have the following format:
#
# id x y z vx vz color
#
# vx and vz are the calculated velocities in x and z direction
# y is calculated form the depth and number of nodes
#
# The output format for the top nodes is
# x y z vx vy vz
#
#

# constants:

FROM_MYRS_TO_YRS = 1e6
FROM_KM_TO_MM = 1e6
FROM_MM_TO_KM = 1e-6
VERSION = "Version 0.4 - 2011.12.16"

def main(argv):
    print VERSION
    # check command line arguments:
    if len(argv) != 2:
        print "wrong number of command line arguments."
        print "usage: {0} input_filename.txt".format(argv[0])
        sys.exit(1)

    inputFileName = argv[1]
    inputFile = None

    # check if file exists and is readable:
    try:
        inputFile = open(inputFileName)
    except IOError as e:
        print "could not open file: {0}".format(inputFileName)
        print e
        sys.exit(1)

    print "reading file '{0}'...".format(inputFileName)

    # read in parameters from input file:
    depthInKM = float(inputFile.readline())
    numberOfNodes = int(inputFile.readline())
    stepYInKM = depthInKM / (numberOfNodes - 1)

    numOfSamples = int(inputFile.readline())
    interpStepX = float(inputFile.readline())

    allFiles = []

    while True:
        line = inputFile.readline().strip()
        if not line: # no more input left in file
            break

        (timeInMyrs, moveFileName) = line.split()
        allFiles.append((float(timeInMyrs), moveFileName))

    inputFile.close()

    # check if there are any input file names given:

    if len(allFiles) == 0:
        print "no imput files specified!"
        sys.exit(1)

    print "depth in km: {0}, number of nodes: {1}, number of input files: {2}".format(depthInKM, numberOfNodes, len(allFiles))
    print "number of samples {0}, x step for interpolation: {1}".format(numOfSamples, interpStepX)

    (lastTimeInMyrs, lastMoveFileName) = allFiles[0]

    lastMoveFile = None
    currentMoveFile = None
    outputFile = None

    xmin = sys.float_info[0]
    xmax = -xmin
    zmin = sys.float_info[0]

    for (currentTimeInMyrs, currentMoveFileName) in allFiles:
        print "processing time: '{0}', file: '{1}' ...".format(currentTimeInMyrs, currentMoveFileName)

        outputFileName = os.path.splitext(currentMoveFileName)[0] # remove file extension if there is any
        topNodesFileName = outputFileName + "_vel_top.dat" # output file name for top nodes
        outputFileName = outputFileName + "_vel_new.dat" # add output extension

        try:
            lastMoveFile = open(lastMoveFileName)
            currentMoveFile = open(currentMoveFileName)
            outputFile = open(outputFileName,"w")
            topNodesFile = open(topNodesFileName,"w")
        except IOError as e:
            print "error while processing files:"
            print lastMoveFileName
            print currentMoveFileName
            print outputFileName
            print topNodesFileName
            print e
            sys.exit(1)

        outputFile.write("# particle id,x,y,z,vx,vy,vz,color; units: km, mm/year\n")
        topNodesFile.write("# x,y,z,vx,vy,vz; units: km, mm/year\n")

        timeDiff = (lastTimeInMyrs - currentTimeInMyrs) * FROM_MYRS_TO_YRS

        lastMoveFilePoints = {}
        currentMoveFilePoints = {}

        for line in lastMoveFile.readlines():
            (x1,_,z1,_,ptid1)=line.split()
            ptid1=int(ptid1); x1=float(x1); z1=float(z1)
            if ptid1 in lastMoveFilePoints:
                print "point IDs not unique!"
                print "file name: {0:s}".format(lastMoveFileName)
                print "old values: id: {0:d}, x: {1:f}, z: {2:f}".format(ptid1, lastMoveFilePoints[ptid1][0], lastMoveFilePoints[ptid1][1])
                print "new values: id: {0:d}, x: {1:f}, z: {2:f}".format(ptid1, x1, z1)
                print
            else:
                lastMoveFilePoints[ptid1] = (x1, z1)

        lastMoveFile.close()

        for line in currentMoveFile.readlines():
            (x1,_,z1,_,ptid1)=line.split()
            ptid1=int(ptid1); x1=float(x1); z1=float(z1)
            if ptid1 in currentMoveFilePoints:
                print "point IDs not unique!"
                print "file name: {0:s}".format(currentMoveFileName)
                print "old values: id: {0:d}, x: {1:f}, z: {2:f}".format(ptid1, currentMoveFilePoints[ptid1][0], currentMoveFilePoints[ptid1][1])
                print "new values: id: {0:d}, x: {1:f}, z: {2:f}".format(ptid1, x1, z1)
                print
            else:
                currentMoveFilePoints[ptid1] = (x1, z1)

        currentMoveFile.close()

        nodesList = []

        for ptid2, (x2, z2) in currentMoveFilePoints.iteritems():
            if ptid2 not in lastMoveFilePoints:
                print "could not find point ID in previous file!"
                print "file1: " + lastMoveFileName
                print "file2: " + currentMoveFileName
                print "id: {0:d}, x: {1:f}, z: {2:f}".format(ptid2, x2, z2)
                print
                continue

            (x1, z1) = lastMoveFilePoints[ptid2]

            # unit of velocity is mm / year

            vx = 0.0
            vz = 0.0

            # input is in km so we have to convert to mmm
            if timeDiff != 0.0:
                vx = FROM_KM_TO_MM * (x2 - x1) / timeDiff
                vz = FROM_KM_TO_MM * (z2 - z1) / timeDiff

            xOut = x2
            zOut = z2

            # unit of positions is meter
            for y in xrange(numberOfNodes):
                yOut = y * stepYInKM

                outputFile.write("{0:20} {1:20} {2:20} {3:20} {4:20} {5:20} {6:20}\n".format(ptid2,xOut,yOut,zOut,vx,0.0,vz))

            if xOut < xmin:
                xmin = xOut
            elif xOut > xmax:
                xmax = xOut
            if zOut < zmin:
                zmin = zOut

            nodesList.append((xOut, zOut, vx, vz))

        # find top nodes:
        valueRange = xmax - xmin

        validNodes = [None for i in xrange(numOfSamples + 1)]

        # take samples: for every range we look for the highest z value
        for node in nodesList:
            index = int(((node[0] - xmin) / valueRange) * numOfSamples)
            if validNodes[index] == None:
                validNodes[index] = node
            elif node[1] > validNodes[index][1]:
                validNodes[index] = node

        print "len validNodes: ", len(validNodes)

        lastNode = None
        outputNodes = []

        # now we need to interpolate between those ranges
        for node in validNodes:
            if node != None:
                (x2, z2, vx2, vz2) = node

                if lastNode != None:
                    (x1, z1, vx1, vz1) = lastNode
                    xInterp = x1
                    zInterp = z1
                    vxInterp = vx1
                    vzInterp = vz1
                    while xInterp < x2:
                        outputNodes.append((xInterp, zInterp, vxInterp, vzInterp))
                        xInterp = xInterp + interpStepX
                        xFactor = (xInterp - x1) / (x2 - x1)
                        zInterp = z1 + (xFactor * (z2 - z1))
                        vxInterp = vx1 + (xFactor * (vx2 - vx1))
                        vzInterp = vz1 + (xFactor * (vz2 - vz1))


                lastNode = (x2, z2, vx2, vz2)

        print "len outputNodes: ", len(outputNodes)

        for (x2, z2, vx2, vz2) in outputNodes:
            topNodesFile.write("{0:20} {1:20} {2:20} {3:20} {4:20} {5:20}\n".format(x2, 0.0, z2, vx2, 0.0, vz2))

        outputFile.close()
        topNodesFile.close()

        lastMoveFileName = currentMoveFileName
        lastTimeInMyrs = currentTimeInMyrs

    print "you need these values for pecube:"
    print "xmin: {0:f}, zmin: {1:f}".format(xmin, zmin)


if __name__ == "__main__":
    main(sys.argv)

