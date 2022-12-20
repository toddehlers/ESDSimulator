#!/usr/bin/env python

# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

# python imports:
import sys
import getopt
import time

# local imports:
import includes.icecascadeConfig
import includes.childConfig
import includes.pecubeConfig


# some globals vars:
PROGRAM_VERSION = "0.4"

def usage():
    print "GridSim V{0:s}".format(PROGRAM_VERSION)
    print "Please specify output format and the filename for the input file."
    print "Allowed values for format: --pecube, --child, --icecascade"
    print
    print "You can also specify the number of CPUs to use locally: -n or --numcpu="
    print
    print "For example:"
    print "{0:s} --pecube -n6 pecube_template/Pecube.in".format(sys.argv[0])
    print
    print "The syntax for each range value is:"
    print "\t boolean values: T:F"
    print "\t integer values: min:max:step, for ex. 10:25:5, step is optional"
    print "\t floating point values: min:max:step, for ex. -8.93:5.67:0.1, step is optional"

def main(argv):

    # default values
    isIceCascade = False
    isChild = False
    isPecube = False
    simCFG = None
    numCPU = 4

    startTime = time.time()

    if len(argv) < 2:
        usage()
        sys.exit(1)

    try:
        opts, args = getopt.getopt(argv, "icpn:", ["icecascade", "child", "pecube", "numcpu="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-i", "--icecascade"):
            isIceCascade = True
        elif opt in ("-c", "--child"):
            isChild = True
        elif opt in ("-p", "--pecube"):
            isPecube = True
        elif opt in ("-n", "--numcpu="):
            numCPU = int(arg)

    inputFile = args[0]

    # TODO: optimize constructor:
    # use base class for parameters and provide a constructor method
    # bc = BaseConfig(inputFile, PROGRAM_VERSION, numCPU)
    # simCFG = bc.createIceConfig()
    # simCFG = bc.createChildConfig()
    # simCFG = bc.createPecubeConfig()

    if isIceCascade:
        simCFG = includes.icecascadeConfig.IceCascadeConfig(inputFile, PROGRAM_VERSION, numCPU)
    elif isChild:
        simCFG = includes.childConfig.ChildConfig(inputFile, PROGRAM_VERSION, numCPU)
    elif isPecube:
        simCFG = includes.pecubeConfig.PecubeConfig(inputFile, PROGRAM_VERSION, numCPU)

    simCFG.loadConfig()
    simCFG.printVariation()
    simCFG.createSimulations()

    endTime = time.time()
    print "time taken: ", (endTime - startTime), " sec."

if __name__ == "__main__":
    main(sys.argv[1:])

