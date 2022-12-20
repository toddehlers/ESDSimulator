import model
import view

class Controller(object):
    def __init__(self):
        self.screenWidth = 640
        self.screenHeight = 640

        self.model = model.Model(self)
        self.view = view.View(self)

        self.modelLoaded = False

    def run(self):
        self.view.run()

    def quit(self):
        pass

    def setROI(self, lowROI, medROI, highROI):
        self.model.setROI(lowROI, medROI, highROI)

    def loadDEM(self, fileName):
        self.model.loadDEM(fileName)
        self.modelLoaded = True

    def loadGDAL(self, fileName):
        self.model.loadGDAL(fileName)
        self.modelLoaded = True

    def loadXML(self, fileName):
        self.model.laodXML(fileName)

    def setDEMProperties(self, ncols, nrows, cellsize, width, height):
        self.view.setDEMProperties(ncols, nrows, cellsize, width, height)

    def isUTM(self):
        return self.view.isUTM()

    def getNODATA_replacement(self):
        return self.view.getNODATA_replacement()

    def getZFactor(self):
        return self.view.getZFactor()

    def getGeoTIFFBand(self):
        return self.view.getGeoTIFFBand()

    def setDisplayData(self, data, minValue, maxValue):
        self.view.setDisplayData(data, minValue, maxValue)

    def updateDisplay(self, currentValue, endValue, message=""):
        self.view.updateDisplay(currentValue, endValue, message)

    def convertROIToModel(self, x1, y1, x2, y2):
        return self.model.convertROIToModel(x1, y1, x2, y2)

    def convertModelToROI(self, x1, y1, w, h):
        return self.model.convertModelToROI(x1, y1, w, h)

    def checkXW(self, x, w):
        return self.model.checkXW(x, w)

    def checkYH(self, y, h):
        return self.model.checkYH(y, h)

    def exportChild(self, outputFileName, childBorder, fixedBorder):
        self.model.exportChild(outputFileName, childBorder, fixedBorder)

    def exportSTL(self, outputFileName, childBorder, fixedBorder):
        self.model.exportSTL(outputFileName, childBorder, fixedBorder)

    def exportEnvi(self, outputFileName, childBorder, fixedBorder):
        self.model.exportEnvi(outputFileName, childBorder, fixedBorder)

    def exportTecplot(self, outputFileName):
        self.model.exportTecplot(outputFileName)

    def exportPecube(self, outputFileName):
        self.model.exportPecube(outputFileName)

    def exportCascade(self, outputFileName):
        self.model.exportCascade(outputFileName)

    def exportCSV(self, outputFileName):
        self.model.exportCSV(outputFileName)

    def getPertubeFactor(self):
        return self.view.getPertubeFactor()

