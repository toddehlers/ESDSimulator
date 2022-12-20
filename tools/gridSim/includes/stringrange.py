# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

import objrange

class StringRange(objrange.ObjRange):
    def __init__(self, name, numOfParamInLine=1, sizeVar=None, sizeVarMod=0):
        objrange.ObjRange.__init__(self, name, numOfParamInLine, sizeVar, sizeVarMod)

    def readLineSpecific(self, values):
        if len(values) != 1:
            raise ValueError(values)

        self.startValue = values[0]
        self.currentValue = self.startValue

    def copy(self):
        return StringRange(self.name, self.numOfParamInLine)

