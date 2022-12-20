# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

import os
import stat
import sys

import boolrange
import stringrange
import intrange
import floatrange

class BaseConfig(object):
    def __init__(self, filename, version, numCPU):
        self.filename = filename

        self.parameters = []
        self.variationParameters = []
        self.stepCounter = 0
        self.outputDir = ""
        self.executable = ""
        self.shellScriptNames = []
        self.inputFileAsParameter = False
        self.currentInputFile = ""
        self.numOfCPU = numCPU
        self.version = version

    def loadConfig(self):
        pass

    def printVariation(self):
        totalSteps = 1
        for parameter in self.parameters:
            if parameter.endValue != None:
                if parameter.step == None:
                    print "variation of parameter '{0:s}' from '{1:s}' to '{2:s}' -> 2 steps".format(parameter.name,
                        str(parameter.startValue),
                        str(parameter.endValue))
                else:
                    print "variation of parameter '{0:s}' from '{1:s}' to '{2:s}' with step '{3:s}' -> {4:d} steps".format(parameter.name,
                        str(parameter.startValue),
                        str(parameter.endValue),
                        str(parameter.step),
                        parameter.numOfSteps)
                totalSteps = totalSteps * parameter.numOfSteps

        print "total number of steps: '{0:d}' (= number of Simulation)".format(totalSteps)

        if totalSteps > 100:
            print "Do you really want to generate more than 100 simulation ? (type 'yes' or 'no')"
            answer = raw_input("---> ")
            if answer != "yes":
                print "You did not enter 'yes', so I will abort now."
                sys.exit(0)

    def createSimulations(self):
        self.stepCounter = 0
        self.variationParameters = []

        print "search variation parameters"

        for parameter in self.parameters:
            if parameter.endValue != None:
                self.variationParameters.append(parameter)

        if not os.access(self.outputDir, os.F_OK):
            print "creating directory: '{0:s}'".format(self.outputDir)
            os.mkdir(self.outputDir)
        else:
            print "using existing directory: '{0:s}'".format(self.outputDir)

        print "numbers of parameters to variate: '{0:d}'".format(len(self.variationParameters))

        self.variateParameters(self.variationParameters)

        masterShellScriptName = "run_all_simulations.sh"
        masterShellScriptFile = open(masterShellScriptName, "w")

        cpuShellScriptFiles = []

        masterShellScriptFile.write("#!/bin/bash\n\n")
        masterShellScriptFile.write("# VERSION: \n\n")
        masterShellScriptFile.write("# This will run all the simulations on this computer locally\n")
        masterShellScriptFile.write("# only {0:d} CPUs are used at the same time\n\n".format(self.numOfCPU))

        for i in range(self.numOfCPU):
            cpuShellScriptFile = "cpu{0:d}.sh".format(i)
            cpuShellScriptFiles.append(cpuShellScriptFile)
            masterShellScriptFile.write("./{0:s} &\n".format(cpuShellScriptFile))

        masterShellScriptFile.close()

        fileInfo = os.stat(masterShellScriptName)
        oldMode = fileInfo.st_mode
        os.chmod(masterShellScriptName, oldMode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        cpuFiles = []

        for cpuShellScriptFile in cpuShellScriptFiles:
            cpuFile = open(cpuShellScriptFile, "w")
            cpuFiles.append(cpuFile)
            cpuFile.write("#!/bin/bash\n")
            cpuFile.write("# This will run the simulations in sequence, not parallel\n")


        fileCounter = 0

        for scriptName in self.shellScriptNames:
            cpuFiles[fileCounter].write("./{0:s} \n".format(scriptName))
            fileCounter = fileCounter + 1
            if fileCounter == len(cpuFiles):
                fileCounter = 0

        for cpuFile in cpuFiles:
            cpuFile.close()

        for cpuShellScriptFile in cpuShellScriptFiles:
            fileInfo = os.stat(cpuShellScriptFile)
            oldMode = fileInfo.st_mode
            os.chmod(cpuShellScriptFile, oldMode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def writeConfigFile(self, outputFile):
        pass

    def writeParameters(self):
        pass

    def writeQsubFile(self):
        qsubName = "{0:s}/{1:04d}/qsub.sh".format(self.outputDir,self.stepCounter)
        qsubFile = open(qsubName, "w")

        qsubFile.write("#!/bin/bash\n\n")
        qsubFile.write("### Options for the cluster simluations:\n")
        qsubFile.write("\n")
        qsubFile.write("### Resouce: one node, one processor per node\n")
        qsubFile.write("#PBS -l nodes=1:ppn=1\n")
        qsubFile.write("\n")
        qsubFile.write("### Make sure bash is default shell\n")
        qsubFile.write("#PBS -S /bin/bash\n")
        qsubFile.write("\n")
        qsubFile.write("### Name of job queue\n")
        qsubFile.write("### - not used now -\n")
        qsubFile.write("\n")
        qsubFile.write("### Name of job\n")
        qsubFile.write("#PBS -N {0:s}_{1:04d}\n".format(self.executable, self.stepCounter))
        qsubFile.write("\n")
        qsubFile.write("### User mail address (change to your needs)\n")
        qsubFile.write("#PBS -M willi.kappler@uni-tuebingen.de\n")
        qsubFile.write("\n")
        qsubFile.write("### Mail options: A = abort, B = begin, E = end\n")
        qsubFile.write("#PBS -m abe\n")
        qsubFile.write("\n")
        qsubFile.write("### Run time of job (aprox)\n")
        qsubFile.write("#PBS -l walltime=120:0:0\n")
        qsubFile.write("\n")
        qsubFile.write("### Go to the correct directory\n")
        qsubFile.write("cd {0:s}/{1:04d}\n".format(self.outputDir,self.stepCounter))
        qsubFile.write("\n")
        qsubFile.write("### Run the executable\n")
        if self.inputFileAsParameter:
            qsubFile.write("./{0:s} {1:s}\n".format(self.executable, self.currentInputFile))
        else:
            qsubFile.write("./{0:s}\n".format(self.executable))
        qsubFile.write("\n")
        qsubFile.write("### Put output into a tar ball\n")
        qsubFile.write("tar cfvz output_{0:s}_{1:04d}.tgz output/\n".format(self.executable, self.stepCounter))
        qsubFile.write("\n")
        qsubFile.write("### Copy archive to file server\n")
        qsubFile.write("### - not used now -\n")
        qsubFile.write("#scp -P6307 output_{0:s}_{1:04d}.tgz result@agassiz.geowissenschaften.uni-tuebingen.de:/data/esd07/cluster/\n".format(self.executable, self.stepCounter))

        qsubFile.close()


        for parameter in self.variationParameters:
            print "parameter: '{0:s}', value: '{1:s}'".format(parameter.name, str(parameter.currentValue))

    def writeShellScripts(self):
        shellScriptName = "{0:s}/{1:04d}/run_this_simulation.sh".format(self.outputDir,self.stepCounter)
        shellScriptFile = open(shellScriptName, "w")

        shellScriptFile.write("#!/bin/bash\n\n")
        shellScriptFile.write("# Change to the coresponding directory\n")
        shellScriptFile.write("cd {0:s}/{1:04d}\n".format(self.outputDir,self.stepCounter))
        shellScriptFile.write("# Run the executable locally on this computer\n")
        shellScriptFile.write("echo 'Running {0:s} in' $(pwd)\n".format(self.executable))
        if self.inputFileAsParameter:
            shellScriptFile.write("./{0:s} {1:s} >& {0:s}_out.txt\n".format(self.executable, self.currentInputFile))
        else:
            shellScriptFile.write("./{0:s} >& {0:s}_out.txt\n".format(self.executable))
        shellScriptFile.write("echo 'Finished {0:s} in' $(pwd)\n".format(self.executable, self.stepCounter))

        shellScriptFile.close()

        fileInfo = os.stat(shellScriptName)
        oldMode = fileInfo.st_mode
        os.chmod(shellScriptName, oldMode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        self.shellScriptNames.append(shellScriptName)

    def variateParameters(self, parameters):
        if len(parameters) == 0:
            self.writeParameters()
            self.writeQsubFile()
            self.writeShellScripts()
            return

        currentParameter = parameters[0]
        recursionList = parameters[1:]

        currentParameter.reset()
        while currentParameter.inRange():
            self.variateParameters(recursionList) # depth first search
            currentParameter.increment()

