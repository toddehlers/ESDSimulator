# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

import objrange

class IntRange(objrange.ObjRange):
    def __init__(self, name, numOfParamInLine=1, sizeVar=None, sizeVarMod=0):
        objrange.ObjRange.__init__(self, name, numOfParamInLine, sizeVar, sizeVarMod)

    def readLineSpecific(self, values):
        numOfValues = len(values)

        if numOfValues == 0 or numOfValues > 3:
            raise ValueError(values)

        self.startValue = int(values[0])

        if numOfValues > 1:
            self.startValue = int(values[0])
            self.endValue = int(values[1])
            if self.startValue >= self.endValue:
                raise ValueError(values)

        if numOfValues == 3:
            self.step = int(values[2])
            if self.step < 0:
                raise ValueError[values]
            self.numOfSteps = ((self.endValue - self.startValue) / self.step) + 1

        self.currentValue = self.startValue

    def copy(self):
        return IntRange(self.name, self.numOfParamInLine)

