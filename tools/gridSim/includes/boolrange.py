# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

import objrange

class BoolRange(objrange.ObjRange):
    def __init__(self, name, numOfParamInLine=1, sizeVar=None, sizeVarMod=0):
        objrange.ObjRange.__init__(self, name, numOfParamInLine, sizeVar, sizeVarMod)

    def readLineSpecific(self, values):
        numOfValues = len(values)
        if numOfValues == 1:
            self.startValue = self.convertValue(values[0])
        elif numOfValues == 2:
            self.startValue = self.convertValue(values[0])
            self.endValue = self.convertValue(values[1])
            if self.startValue == self.endValue:
                raise ValueError(values)
        else:
            raise ValueError(values)

        self.currentValue = self.startValue

    def convertValue(self, value):
        if value == "T":
            return True
        elif value == "F":
            return False
        else:
            raise ValueError(value)

    def valueAsString(self):
        if self.currentValue:
            return "T"
        else:
            return "F"

    def copy(self):
        return BoolRange(self.name, self.numOfParamInLine)

