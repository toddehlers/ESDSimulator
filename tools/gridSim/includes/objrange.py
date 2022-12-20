# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

import math

class ObjRange(object):
    def __init__(self, name, numOfParamInLine, sizeVar, sizeVarMod):
        self.name = name
        self.numOfParamInLine = numOfParamInLine
        self.sizeVar = sizeVar
        self.sizeVarMod = sizeVarMod
        self.sizeVarValue = None
        self.startValue = None
        self.endValue = None
        self.currentValue = None
        self.step = None
        self.numOfSteps = None
        self.comments = []

    def readLine(self, line):
        self.readLineSpecific(line.split(":"))
#        print "read line for '{0:s}', '{1:s}', '{2:s}', '{3:s}', '{4:s}', '{5:s}'".format(self.name,
#                line,
#                self.__class__.__name__,
#                str(self.currentValue),
#                str(self.endValue),
#                str(self.step))

        if self.step == None:
            if self.endValue != None:
                self.numOfSteps = 2
            else:
                self.numOfSteps = 1

    def readLineSpecific(self, values):
        pass

    def reset(self):
        self.currentValue = self.startValue

    def setEndValue(self):
        self.currentValue = self.endValue

    def inRange(self):
        if self.step == None:
            return self.currentValue != None
        else:
            return self.currentValue <= self.endValue

    def increment(self):
        if self.step == None:
            if self.currentValue == self.startValue:
                self.currentValue = self.endValue
            elif self.currentValue == self.endValue:
                self.currentValue = None
        else:
            self.currentValue = self.currentValue + self.step

    def valueAsString(self):
        return str(self.currentValue)

    def copy(self):
        pass

