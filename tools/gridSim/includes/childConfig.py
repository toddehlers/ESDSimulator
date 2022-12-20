# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

# system import
import os
import sys

# local import
import baseConfig
import boolrange
import stringrange
import intrange
import floatrange

class ChildConfig(baseConfig.BaseConfig):
    def __init__(self, filename, version, numCPU):
        baseConfig.BaseConfig.__init__(self, filename, version, numCPU)

        self.outputDir = "child"
        self.executable = "child"
        self.inputFileAsParameter = True

    def loadConfig(self):
        inFile = open(self.filename, "r")
        comments = []

        while True:
            line = inFile.readline()
            if line == "":
                break

            line = line.strip()
            if line.startswith("#") or len(line) == 0:
                comments.append(line)
            else:
                values = line.split(":", 1)
                paramName = values[0].strip()

                if len(values) > 1:
                    paramDescription = values[1].strip()
                else:
                    paramDescription = ""

                nextLine = inFile.readline().strip()
                values = nextLine.split(":")

                # try to figure out what type the input is

                newParameter = None

                try:
                    maybeInt = int(values[0])
                    newParameter = intrange.IntRange(paramName)

                except ValueError:
                    # not an integer, maybe a float ?
                    try:
                        maybeFloat = float(values[0])
                        newParameter = floatrange.FloatRange(paramName)

                    except ValueError:
                        # ok, neither int nor float, so it must be a string!
                        newParameter = stringrange.StringRange(paramName)

                newParameter.comments = comments
                newParameter.description = paramDescription

#                print "name: '{0:s}', desc: '{1:s}', nextLine: '{2:s}'".format(newParameter.name, newParameter.description, nextLine)

                newParameter.readLine(nextLine)

                self.parameters.append(newParameter)

                comments = []

        inFile.close()

    def writeConfigFile(self, outputFile):
        for parameter in self.parameters:
            for comment in parameter.comments:
                outputFile.write(comment + "\n")
            outputFile.write("{0:s}: {1:s}\n".format(parameter.name, parameter.description))
            outputFile.write(parameter.valueAsString() + "\n")

    def writeParameters(self):
        self.stepCounter = self.stepCounter + 1
        print "stepCounter: ", self.stepCounter

        currentDirectory = "{0:s}/{1:04d}".format(self.outputDir, self.stepCounter)

        try:
            os.mkdir(currentDirectory)
        except OSError:
            pass # using existing directory

        os.symlink("../../child_template/child", "{0:s}/child".format(currentDirectory))
        os.symlink("../../child_template/input", "{0:s}/input".format(currentDirectory))

        self.currentInputFile = os.path.basename(self.filename)
        paramOutFile = open("{0:s}/{1:s}".format(currentDirectory, self.currentInputFile), "w")

        self.writeConfigFile(paramOutFile)

        paramOutFile.close()

