import math
import sys
import random
import time
import os.path
import struct

from osgeo import gdal
from osgeo import gdalconst

# import matplotlib.delaunay as triang

# WK 2011.07.15: for debugging:
#import guppy
#import objgraph
#import gc

import main

EARTH_RADIUS_IN_M = 6371000.785
EARTH_CIRCUMFENCE_IN_M = EARTH_RADIUS_IN_M * 2 * math.pi
ONE_DEGREE_IN_M = EARTH_CIRCUMFENCE_IN_M / 360.0 # unit: meter / degree

class Point(object):
    def __init__(self, px, py, pz, b):
        self.x = px
        self.y = py
        self.z = pz
        self.border = b

class ZoneROI(object):
    def __init__(self, outerROI, innerROI):
        self.outerROI = outerROI
        self.innerROI = innerROI

        self.step = (self.outerROI.getStep() + self.innerROI.getStep()) / 2

        (x1, y1, w1, h1) = self.outerROI.getModelROI()
        (x2, y2, w2, h2) = self.innerROI.getModelROI()

        outerFactor = 0.2
        innerFactor = 0.8

        self.x = (outerFactor * x1) + (innerFactor * x2)
        self.y = (outerFactor * y1) + (innerFactor * y2)
        self.w = (outerFactor * w1) + (innerFactor * w2)
        self.h = (outerFactor * h1) + (innerFactor * h2)

    def getStep(self):
        return self.step

    def getModelROI(self):
        return (self.x, self.y, self.w, self.h)

class Model(object):
    def __init__(self, controller):
        self.controller = controller
        self.DEMFileName = None
        self.XMLFileName = None
        self.lowROI = None
        self.medROI = None
        self.highROI = None
        self.GDALUsed = False

    def loadDEM(self, fileName):
        self.DEMFileName = fileName

        DEMFile = open(self.DEMFileName)

        self.ncols, self.nrows, self.cellsize, self.NODATA_value = self.loadHeader(DEMFile)

        if self.controller.isUTM():
            self.UTMFactor = 0.001
        else:
            self.UTMFactor = ONE_DEGREE_IN_M / 1000.0

        self.width = self.ncols * self.UTMFactor * self.cellsize
        self.height = self.nrows * self.UTMFactor * self.cellsize
        self.cellsizeInKilometer = self.cellsize * self.UTMFactor

        self.controller.setDEMProperties(self.ncols, self.nrows, self.cellsize, self.width, self.height)

        scaleFactor = max(self.ncols / 640, self.nrows / 640) + 1
        print "scaleFactor: ", scaleFactor

        rowCounter = 0

        rowData = []

        self.NODATA_replacement = self.controller.getNODATA_replacement()
        self.zFactor = self.controller.getZFactor()

        minValue = sys.float_info[0]
        maxValue = sys.float_info[3]

        for i in xrange(self.nrows):
            line = DEMFile.readline()
            rowCounter = rowCounter + 1
            if rowCounter >= scaleFactor:
                rowCounter = 0
                colCounter = 0
                colData = []
                lineValues = line.split()
                self.controller.updateDisplay(i, self.nrows, "loadDEM")
                for j in xrange(self.ncols):
                    colCounter = colCounter + 1
                    if colCounter >= scaleFactor:
                        colCounter = 0

                        value = self.convertValue(lineValues[j])

                        if value > maxValue:
                            maxValue = value
                        elif value < minValue:
                            minValue = value
                        colData.append(value)
                        if len(colData) >= 640:
                            break
                rowData.append(colData)
                if len(rowData) >= 640:
                    break

        DEMFile.close()

        self.GFXRows = len(rowData)
        self.GFXCols = len(rowData[0])

        print "rowData: ", self.GFXRows
        print "colData: ", self.GFXCols
        print "min, max: ", minValue, maxValue

        self.controller.setDisplayData(rowData, minValue, maxValue)
        self.GDALUsed = False

    def loadGDAL(self, fileName):
        self.DEMFileName = fileName

        if self.controller.isUTM():
            self.UTMFactor = 0.001
        else:
            self.UTMFactor = ONE_DEGREE_IN_M / 1000.0

        dataset = gdal.Open(self.DEMFileName, gdalconst.GA_ReadOnly)

        self.ncols = dataset.RasterXSize
        self.nrows = dataset.RasterYSize

        geotransform = dataset.GetGeoTransform()

        self.cellsize = abs(geotransform[1])

        self.width = self.ncols * self.UTMFactor * self.cellsize
        self.height = self.nrows * self.UTMFactor * self.cellsize
        self.cellsizeInKilometer = self.cellsize * self.UTMFactor

        self.controller.setDEMProperties(self.ncols, self.nrows, self.cellsize, self.width, self.height)

        scaleFactor = max(self.ncols / 640, self.nrows / 640) + 1
        print "scaleFactor: ", scaleFactor

        rowCounter = 0

        rowData = []

        self.NODATA_replacement = self.controller.getNODATA_replacement()
        self.zFactor = self.controller.getZFactor()

        minValue = sys.float_info[0]
        maxValue = sys.float_info[3]

        rasterCount = dataset.RasterCount
        currentBand = self.controller.getGeoTIFFBand()

        if not currentBand in xrange(1, rasterCount + 1):
            currentBand = 1

        dataBand = dataset.GetRasterBand(currentBand)

        for i in xrange(self.nrows):
            rowCounter = rowCounter + 1
            if rowCounter >= scaleFactor:
                rowCounter = 0
                colCounter = 0
                colData = []
                self.controller.updateDisplay(i, self.nrows, "loadGDAL")
                for j in xrange(self.ncols):
                    colCounter = colCounter + 1
                    if colCounter >= scaleFactor:
                        colCounter = 0

                        # TODO: use scan line from ReadRaster

                        data = dataBand.ReadRaster(j, i, 1, 1, 1, 1, gdalconst.GDT_Float32)
                        data = struct.unpack('f', data)
                        value = data[0]

                        if value > maxValue:
                            maxValue = value
                        elif value < minValue:
                            minValue = value
                        colData.append(value * self.zFactor)
                        if len(colData) >= 640:
                            break
                rowData.append(colData)
                if len(rowData) >= 640:
                    break

        self.GFXRows = len(rowData)
        self.GFXCols = len(rowData[0])

        print "rowData: ", self.GFXRows
        print "colData: ", self.GFXCols
        print "min, max: ", minValue, maxValue

        self.controller.setDisplayData(rowData, minValue, maxValue)
        self.GDALUsed = True

    def loadHeader(self, DEMFile):
        # header of input file looks like this:
        #
        # ncols         6893            #
        # nrows         6429            #
        # xllcorner     0.32568359375   #
        # yllcorner     41.884765625    #
        # cellsize      0.00208333333   #
        # NODATA_value  -9999           #
        #

        ncols = self.readInt(DEMFile)
        nrows = self.readInt(DEMFile)
        self.readFloat(DEMFile)
        self.readFloat(DEMFile)
        cellsize = self.readFloat(DEMFile)
        NODATA_value = self.readFloat(DEMFile)

        return (ncols, nrows, cellsize, NODATA_value)

    def readInt(self, inFile):
        line = inFile.readline()
        content = line.split()
        if len(content) < 2:
            raise ValueError(line, "invalid data")
        try:
            value = int(content[1])
            return value
        except ValueError:
            raise ValueError(line, "not an int")

    def readFloat(self, inFile):
        line = inFile.readline()
        content = line.split()
        try:
            value = float(content[1])
            return value
        except ValueError:
            raise ValueError(line, "not a float")

    def setROI(self, lowROI, medROI, highROI):
        self.lowROI = lowROI
        self.medROI = medROI
        self.highROI = highROI

    def loadXML(self, fileName):
        self.XMLFileName = fileName

    def checkXW(self, x, w):
        if x < 0.0:
            x = 0.0

        if x >= self.width:
            if w > self.width:
                x, w = 0.0, self.width
            else:
                x = self.width - w
        else:
            if x+w > self.width:
                w = self.width - x

        return (x, w)

    def checkYH(self, y, h):
        if y < 0.0:
            y = 0.0

        if y >= self.height:
            if h > self.height:
                y, h = 0.0, self.height
            else:
                y = self.height - h
        else:
            if y+h > self.height:
                h = self.height - y

        return (y, h)

    def ROIXToModel(self, value):
        modelX = value * self.width / self.GFXCols
        if modelX > self.width:
            return self.width
        return modelX

    def modelXToROI(self, value):
        ROIX = value * self.GFXCols / self.width
        if ROIX > self.controller.screenWidth:
            ROIX = self.controller.screenWidth
        return int(ROIX)

    def ROIYToModel(self, value):
        modelY = value * self.height / self.GFXRows
        if modelY > self.height:
            return self.height
        return modelY

    def modelYToROI(self, value):
        ROIY = value * self.GFXRows / self.height
        if ROIY > self.controller.screenHeight:
            ROIY = self.controller.screenHeight
        return int(ROIY)

    def convertROIToModel(self, x1, y1, x2, y2):
        mx1 = self.ROIXToModel(x1)
        my1 = self.ROIYToModel(y1)
        mx2 = self.ROIXToModel(x2)
        my2 = self.ROIYToModel(y2)
        return (mx1, my1, mx2 - mx1, my2 - my1)

    def convertModelToROI(self, x1, y1, w, h):
        rx1 = self.modelXToROI(x1)
        ry1 = self.modelYToROI(y1)
        rx2 = self.modelXToROI(x1 + w)
        ry2 = self.modelYToROI(y1 + h)
        return (rx1, ry1, rx2, ry2)

    def linearInterpolation(self, a, b, c, d, e):
        return ((c-a)*d/(c-b)) + ((a-b)*e/(c-b))

    def bilinearInterpolation(self, p1, p2, p3, p4, nx, ny):
        r1 = self.linearInterpolation(nx,p1.x,p2.x,p3.z,p4.z)
        r2 = self.linearInterpolation(nx,p1.x,p2.x,p1.z,p2.z)
        return self.linearInterpolation(ny,p3.y,p1.y,r1,r2)

    def exportChild(self, outputFilename, childBorder, fixedBorder):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportChild")

        if self.lowROI.exportROI():
            self.setChildBorder(lowData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.medROI.exportROI():
            self.setChildBorder(medData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.highROI.exportROI():
            self.setChildBorder(highData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        else:
            return

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportChild")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        childFile = open(outputFilename, "w")
        childFile.write("{0:d}\n".format(len(allData)))

        for point in allData:
            childFile.write("{0:f} {1:f} {2:f} {3:d}\n".format(point.x * 1000.0, (self.height - point.y) * 1000.0, point.z, point.border))

        childFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportSTL(self, outputFilename, childBorder, fixedBorder):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportSTL")

        if self.lowROI.exportROI():
            self.setChildBorder(lowData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.medROI.exportROI():
            self.setChildBorder(medData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.highROI.exportROI():
            self.setChildBorder(highData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        else:
            return

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportSTL")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        STLFile = open(outputFilename, "w")

        # write STL header

        STLFile.write("solid surface\n")

        coordinateX = []
        coordinateY = []

        for point in allData:
            coordinateX.append(point.x * 1000.0)
            coordinateY.append(point.y * 1000.0)

        triangles = []

# 2013.07.11, WK: does not work an gagssiz, ELF32 liblapack needs to be removed
#        cens,edg,triangles,neig = triang.delaunay(coordinateX, coordinateY)

        for tri in triangles:
            tri1 = allData[tri[0]]
            tri2 = allData[tri[1]]
            tri3 = allData[tri[2]]

            STLFile.write("\tfacet normal 0.0 0.0 0.0\n")
            STLFile.write("\t\touter loop\n")
            STLFile.write("\t\t\tvertex {0:f} {1:f} {2:f}\n".format(tri1.x * 1000.0, tri1.y * 1000.0, tri1.z))
            STLFile.write("\t\t\tvertex {0:f} {1:f} {2:f}\n".format(tri2.x * 1000.0, tri2.y * 1000.0, tri2.z))
            STLFile.write("\t\t\tvertex {0:f} {1:f} {2:f}\n".format(tri3.x * 1000.0, tri3.y * 1000.0, tri3.z))
            STLFile.write("\t\tendloop\n")
            STLFile.write("\tendfacet\n")


        STLFile.write("endsolid\n")

        STLFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportEnvi(self, outputFilename, childBorder, fixedBorder):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportSTL")

        if self.lowROI.exportROI():
            self.setChildBorder(lowData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.medROI.exportROI():
            self.setChildBorder(medData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        elif self.highROI.exportROI():
            self.setChildBorder(highData, childBorder[0], childBorder[1], childBorder[2], childBorder[3], fixedBorder)
        else:
            return

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportSTL")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        # Insert GDAL Code here




        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportTecplot(self, outputFilename):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportTecplot")

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportTecplot")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        tecplotFile = open(outputFilename, "w")
        tecplotFile.write('TITLE = "convertArc export, version: {0:s}"\n'.format(main.VERSION))
        tecplotFile.write('VARIABLES = "x" "y" "z"\n')
        tecplotFile.write('ZONE T = "{0:s}"\n'.format(self.DEMFileName))
        tecplotFile.write('i = {0:d}, j = 1, k = 1, f = point\n'.format(len(allData)))

# Zone I=78, J=78, K=1, F=POINT

        print "nodes: {0:d}".format(len(allData))
#        print "point: {0:s}".format(str(allData[0]))

        for point in allData:
            tecplotFile.write("{0:f} {1:f} {2:f}\n".format(point.x, self.height - point.y, point.z))

        tecplotFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportCSV(self, outputFilename):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportCSV")

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportCSV")

        csvFile = open(outputFilename, "w")
        csvFile.write("x, y, z\n")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        print "nodes: {0:d}".format(len(allData))

        for point in allData:
            csvFile.write("{0:f}, {1:f}, {2:f}\n".format(point.x, self.height - point.y, point.z))

        csvFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportCascade(self, outputFilename):
        startTime = time.time()
        self.setModelConditions()

        (lowData, medData, highData, lowZoneData, medZoneData) = self.retrieveAllModelData("exportCascade")

        (lowData, medData, highData, lowZoneData, medZoneData) = self.perturbeAllData(lowData, medData, highData, lowZoneData, medZoneData, "exportCascade")

        allData = lowData + medData + highData + lowZoneData + medZoneData

        cascadeFile = open(outputFilename, "w")
        cascadeFile.write('TITLE = "convertArc export, version: {0:s}"\n'.format(main.VERSION))
        cascadeFile.write('VARIABLES = "x" "y" "z" "node" "Precipitation" "Fluvial Erosion Rate" "Diffusion Erosion Rate" "Landslide Erosion Rate" "Total Erosion Rate" "Catchment Color" "Catchment Number" "Glacial Erosion Rate" "Ice Thickness" "Mass Balance" "Total Topography" "Sliding Velocity" "Sediment flux" "Rock/Alluvium Contact" "Isostatic deflection" "Slope" "Totalflexiso" "Constriction"\n')
        cascadeFile.write('ZONE T = "{0:s}"\n'.format(self.DEMFileName))
        cascadeFile.write('i = {0:d}, j = 1, k = 1, f = point\n'.format(len(allData)))

        print "nodes: {0:d}".format(len(allData))

        for point in allData:
            cascadeFile.write("{0:f} {1:f} {2:f} 0 0.0 0.0 0.0 0.0 0.0 0 0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n".format(point.x, self.height - point.y, point.z))

        cascadeFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def exportPecube(self, outputFilename):
        startTime = time.time()
        self.setModelConditions()

        # TODO: export med or high data to pecube
        if self.lowROI.exportROI():
            lowData, nodes = self.retrieveModelData(self.yInROI(self.lowROI),
                                                    self.xInROI(self.lowROI),
                                                    0, "exportPecube retrieve lowData", self.lowROI.getStep())
        else:
            return

        print "nodes: {0:d}".format(nodes)

        print "size of Point in data: ", sys.getsizeof(lowData[0][0])
        print "size of line in data: ", sys.getsizeof(lowData[0])
        print "size of data: ", sys.getsizeof(lowData)


# WK 2011.07.15: for debugging:
#        print "heap: ", guppy.hpy().heap()
#        objgraph.show_refs(lowData, filename='lowData.png')

        numOfPointsX = len(lowData[0])
        numOfPointsY = len(lowData)

        pecubeFile = open(outputFilename, "w")

        for line in reversed(lowData):
            for point in line:
                pecubeFile.write("{0:f}\n".format(point.z))

        pecubeFile.close()

        (head, tail) = os.path.splitext(outputFilename)

        pecubeHeaderFile = open(head + "_header" + tail, "w")
        pecubeHeaderFile.write("number of points in x: {0:d} (Input 6)\n".format(numOfPointsX))
        pecubeHeaderFile.write("number of points in y: {0:d} (Input 6)\n".format(numOfPointsY))
        pecubeHeaderFile.write("spacing in x: {0:f} (Input 7)\n".format(self.cellsizeInKilometer * 1000.0))
        pecubeHeaderFile.write("spacing in y: {0:f} (Input 7)\n".format(self.cellsizeInKilometer * 1000.0))
        pecubeHeaderFile.close()

        endTime = time.time()

        print "time taken: {0:f} sec.".format(endTime - startTime)

    def setModelConditions(self):
        self.zFactor = self.controller.getZFactor()
        self.NODATA_replacement = self.controller.getNODATA_replacement()
        self.perturbeFactor = self.controller.getPertubeFactor()
        print "perturbe factor: {0:f}".format(self.perturbeFactor)

    def xInROI(self, outerROI, innerROI=None):
        outerROIX1, outerROIY1, width, height = outerROI.getModelROI()
        outerROIX2 = outerROIX1 + width
        outerROIY2 = outerROIY1 + height

        stepSize = outerROI.getStep()

#        print "xInROI, outerROI old values: {0:f}, {1:f}, {2:f}, {3:f}".format(outerROIX1, outerROIY1, outerROIX2, outerROIY2)

        outerROIX1 = outerROIX1 - (stepSize * self.cellsizeInKilometer)
        outerROIY1 = outerROIY1 - (stepSize * self.cellsizeInKilometer)
        outerROIX2 = outerROIX2 + (stepSize * self.cellsizeInKilometer)
        outerROIY2 = outerROIY2 + (stepSize * self.cellsizeInKilometer)

#        print "xInROI, outerROI new values: {0:f}, {1:f}, {2:f}, {3:f}\n".format(outerROIX1, outerROIY1, outerROIX2, outerROIY2)

        if innerROI == None:
            def oneROI(currentXInKM, currentYInKM):
                return currentXInKM > outerROIX1 and currentXInKM < outerROIX2

            return oneROI

        innerROIX1, innerROIY1, width, height = innerROI.getModelROI()
        innerROIX2 = innerROIX1 + width
        innerROIY2 = innerROIY1 + height

#        print "xInROI, innerROI old values: {0:f}, {1:f}, {2:f}, {3:f}".format(innerROIX1, innerROIY1, innerROIX2, innerROIY2)

        innerROIX1 = innerROIX1 + (stepSize * self.cellsizeInKilometer)
        innerROIY1 = innerROIY1 + (stepSize * self.cellsizeInKilometer)
        innerROIX2 = innerROIX2 - (stepSize * self.cellsizeInKilometer)
        innerROIY2 = innerROIY2 - (stepSize * self.cellsizeInKilometer)

#        print "xInROI, innerROI new values: {0:f}, {1:f}, {2:f}, {3:f}\n".format(innerROIX1, innerROIY1, innerROIX2, innerROIY2)

        def twoROIS(currentXInKM, currentYInKM):
            if currentXInKM > outerROIX1 and currentXInKM < outerROIX2:
                if currentXInKM < innerROIX1 or currentXInKM > innerROIX2:
                    return True
                elif currentYInKM < innerROIY1 or currentYInKM > innerROIY2:
                    return True
            else:
                return False

        return twoROIS

    def yInROI(self, outerROI, innerROI = None):
        outerROIX1, outerROIY1, width, height = outerROI.getModelROI()
        outerROIX2 = outerROIX1 + width
        outerROIY2 = outerROIY1 + height

        stepSize = outerROI.getStep()

#        print "yInROI, outerROI old values: {0:f}, {1:f}, {2:f}, {3:f}".format(outerROIX1, outerROIY1, outerROIX2, outerROIY2)

        outerROIX1 = outerROIX1 - (stepSize * self.cellsizeInKilometer)
        outerROIY1 = outerROIY1 - (stepSize * self.cellsizeInKilometer)
        outerROIX2 = outerROIX2 + (stepSize * self.cellsizeInKilometer)
        outerROIY2 = outerROIY2 + (stepSize * self.cellsizeInKilometer)

#        print "yInROI, outerROI new values: {0:f}, {1:f}, {2:f}, {3:f}\n".format(outerROIX1, outerROIY1, outerROIX2, outerROIY2)

        if innerROI == None:
            def oneROI(currentXInKM, currentYInKM):
                return currentYInKM > outerROIY1 and currentYInKM < outerROIY2

            return oneROI

        innerROIX1, innerROIY1, width, height = innerROI.getModelROI()
        innerROIX2 = innerROIX1 + width
        innerROIY2 = innerROIY1 + height

#        print "yInROI, innerROI old values: {0:f}, {1:f}, {2:f}, {3:f}".format(innerROIX1, innerROIY1, innerROIX2, innerROIY2)

        innerROIX1 = innerROIX1 + (stepSize * self.cellsizeInKilometer)
        innerROIY1 = innerROIY1 + (stepSize * self.cellsizeInKilometer)
        innerROIX2 = innerROIX2 - (stepSize * self.cellsizeInKilometer)
        innerROIY2 = innerROIY2 - (stepSize * self.cellsizeInKilometer)

#        print "yInROI, innerROI new values: {0:f}, {1:f}, {2:f}, {3:f}\n".format(innerROIX1, innerROIY1, innerROIX2, innerROIY2)

        def twoROIS(currentXInKM, currentYInKM):
            if currentYInKM > outerROIY1 and currentYInKM < outerROIY2:
                if currentYInKM < innerROIY1 or currentYInKM > innerROIY2:
                    return True
                elif currentXInKM < innerROIX1 or currentXInKM > innerROIX2:
                    return True
            else:
                return False

        return twoROIS

    def convertValue(self, value):
        value = float(value)
        if value == self.NODATA_value:
            return self.NODATA_replacement * self.zFactor
        else:
            return value * self.zFactor

    def retrieveModelData(self, yInROI, xInROI, nodes, GUIMessage, stepSize):
        if self.GDALUsed:
            return self.retrieveModelDataGDAL(yInROI, xInROI, nodes, GUIMessage, stepSize)
        else:
            return self.retrieveModelDataDEM(yInROI, xInROI, nodes, GUIMessage, stepSize)

    def retrieveModelDataDEM(self, yInROI, xInROI, nodes, GUIMessage, stepSize):
        DEMFile = open(self.DEMFileName)

        self.loadHeader(DEMFile)

        xInKM = 0.0
        yInKM = 0.0
        currentXStep = 0
        currentYStep = 0
        data = []

        xStepSize = stepSize
        yStepSize = stepSize

        for row in xrange(self.nrows):
            line = DEMFile.readline()

            currentYStep = currentYStep + 1
            if currentYStep >= yStepSize:
                currentYStep = 0
                if yInROI(xInKM, yInKM):
                    xInKM = 0.0
                    currentXStep = 0
                    line = line.split()
                    lineValues = []
                    self.controller.updateDisplay(row, self.nrows, GUIMessage)

                    for column in xrange(self.ncols):
                        currentXStep = currentXStep + 1
                        if currentXStep >= xStepSize:
                            currentXStep = 0
                            if xInROI(xInKM, yInKM):
                                lineValues.append(Point(xInKM, yInKM, self.convertValue(line[column]), 0))
                                nodes = nodes + 1

                        xInKM = xInKM + self.cellsizeInKilometer

                    data.append(lineValues)

            yInKM = yInKM + self.cellsizeInKilometer

        DEMFile.close()

        return (data, nodes)

    def retrieveModelDataGDAL(self, yInROI, xInROI, nodes, GUIMessage, stepSize):

        xInKM = 0.0
        yInKM = 0.0
        currentXStep = 0
        currentYStep = 0
        data = []

        xStepSize = stepSize
        yStepSize = stepSize

        dataset = gdal.Open(self.DEMFileName, gdalconst.GA_ReadOnly)
        rasterCount = dataset.RasterCount

        currentBand = self.controller.getGeoTIFFBand()

        if not currentBand in xrange(1, rasterCount + 1):
            currentBand = 1

        dataBand = dataset.GetRasterBand(currentBand)

        for row in xrange(self.nrows):

            currentYStep = currentYStep + 1
            if currentYStep >= yStepSize:
                currentYStep = 0
                if yInROI(xInKM, yInKM):
                    xInKM = 0.0
                    currentXStep = 0
                    lineValues = []
                    self.controller.updateDisplay(row, self.nrows, GUIMessage)

                    for column in xrange(self.ncols):
                        currentXStep = currentXStep + 1
                        if currentXStep >= xStepSize:
                            currentXStep = 0
                            if xInROI(xInKM, yInKM):

                                # TODO: make it more efficient and use scan line from ReadRaster

                                dataB = dataBand.ReadRaster(column, row, 1, 1, 1, 1, gdalconst.GDT_Float32)
                                dataB = struct.unpack('f', dataB)

                                lineValues.append(Point(xInKM, yInKM, dataB[0], 0))
                                nodes = nodes + 1

                        xInKM = xInKM + self.cellsizeInKilometer

                    data.append(lineValues)

            yInKM = yInKM + self.cellsizeInKilometer

        return (data, nodes)

    def retrieveAllModelData(self, GUIMessage):
        lowData = []
        medData = []
        highData = []

        lowZoneData = []
        medZoneData = []

        nodeCounter = 0

        if self.lowROI.exportROI():
            if self.medROI.exportROI():

                lowZoneROI = ZoneROI(self.lowROI, self.medROI)

                lowData, nodeCounter = self.retrieveModelData(self.yInROI(self.lowROI, lowZoneROI),
                                                             self.xInROI(self.lowROI, lowZoneROI),
                                                             0, GUIMessage + " retrieve lowData", self.lowROI.getStep())

                lowZoneData, nodeCounter = self.retrieveModelData(self.yInROI(lowZoneROI, self.medROI),
                                                             self.xInROI(lowZoneROI, self.medROI),
                                                             0, GUIMessage + " retrieve lowZoneData", lowZoneROI.getStep())

#                self.medROI.setCoordinatesToGrid(self.lowStep * self.cellsizeInKilometer, True)
            else:
                if self.highROI.exportROI():

                    lowZoneROI = ZoneROI(self.lowROI, self.highROI)

                    lowData, nodeCounter = self.retrieveModelData(self.yInROI(self.lowROI, lowZoneROI),
                                                             self.xInROI(self.lowROI, lowZoneROI),
                                                             0, GUIMessage + " retrieve lowData", self.lowROI.getStep())

                    lowZoneData, nodeCounter = self.retrieveModelData(self.yInROI(lowZoneROI, self.highROI),
                                                             self.xInROI(lowZoneROI, self.highROI),
                                                             0, GUIMessage + " retrieve lowZoneData", lowZoneROI.getStep())
#                    self.highROI.setCoordinatesToGrid(self.lowStep * self.cellsizeInKilometer, True)
                else:
                    lowData, nodeCounter = self.retrieveModelData(self.yInROI(self.lowROI),
                                                                 self.xInROI(self.lowROI),
                                                                 0, GUIMessage + " retrieve lowData", self.lowROI.getStep())

        if self.medROI.exportROI():
            if self.highROI.exportROI():

                medZoneROI = ZoneROI(self.medROI, self.highROI)

                medData, nodeCounter = self.retrieveModelData(self.yInROI(self.medROI, medZoneROI),
                                                             self.xInROI(self.medROI, medZoneROI),
                                                             nodeCounter, GUIMessage + " retrieve medData", self.medROI.getStep())

                medZoneData, nodeCounter = self.retrieveModelData(self.yInROI(medZoneROI, self.highROI),
                                                             self.xInROI(medZoneROI, self.highROI),
                                                             nodeCounter, GUIMessage + " retrieve medZoneData", medZoneROI.getStep())
#                self.highROI.setCoordinatesToGrid(self.medStep * self.cellsizeInKilometer, True)
            else:
                medData, nodeCounter = self.retrieveModelData(self.yInROI(self.medROI),
                                                             self.xInROI(self.medROI),
                                                             nodeCounter, GUIMessage + " retrieve medData", self.medROI.getStep())

        if self.highROI.exportROI():
            highData, nodeCounter = self.retrieveModelData(self.yInROI(self.highROI),
                                                          self.xInROI(self.highROI),
                                                          nodeCounter, GUIMessage + " retrieve highData", self.highROI.getStep())


        print "total nodes: {0:d}".format(nodeCounter)

        return (lowData, medData, highData, lowZoneData, medZoneData)


    def setChildBorder(self, data, top, bottom, left, right, fixedBorder):
        # change values for the first line
        for col in xrange(len(data[0])):
            data[0][col].border = top
            if fixedBorder != None:
                data[0][col].z = fixedBorder

        # change values for the last line
        for col in xrange(len(data[-1])):
            data[-1][col].border = bottom
            if fixedBorder != None:
                data[-1][col].z = fixedBorder

        # change values for the left side
        for row in xrange(len(data)):
            data[row][0].border = left
            if fixedBorder != None:
                data[row][0].z = fixedBorder

        # change values for the right side
        for row in xrange(len(data)):
            data[row][-1].border = right
            if fixedBorder != None:
                data[row][-1].z = fixedBorder

    def perturbeData(self, data, reallyPertube, GUIMessage):
        result = []

        if len(data) < 3:
            print "perturbe data len: {0:d} - {1:s}".format(len(data), GUIMessage)
            return [point for line in data for point in line]

        if len(data[0]) < 3:
            print "perturbe data[0] len: {0:d} - {1:s}".format(len(data[0]), GUIMessage)
            return [point for line in data for point in line]

        pointDistance = data[0][1].x - data[0][0].x
        tolerance = pointDistance / 2.0
        numOfRows = len(data)

        for row in xrange(1, numOfRows - 1):
            self.controller.updateDisplay(row, numOfRows, GUIMessage)
            for col in xrange(1, len(data[row]) - 1):
                currentPoint = data[row][col]
                # find matching points in previous row
                p1, p2 = self.findBiPoints(data[row-1], currentPoint.y - pointDistance,
                                           currentPoint.x - pointDistance,
                                           currentPoint.x + pointDistance,
                                           tolerance)
                # find points in next row
                p3, p4 = self.findBiPoints(data[row+1], currentPoint.y + pointDistance,
                                           currentPoint.x - pointDistance,
                                           currentPoint.x + pointDistance,
                                           tolerance)
                if (p1 != None) and (p2 != None) and (p3 != None) and (p4 != None):
                    if reallyPertube:
                        newX = currentPoint.x + ((random.random() - 0.5) * pointDistance * self.perturbeFactor)
                        newY = currentPoint.y + ((random.random() - 0.5) * pointDistance * self.perturbeFactor)
                        newZ = self.bilinearInterpolation(p1, p2, p3, p4, newX, newY)
                        result.append(Point(newX, newY, newZ, currentPoint.border))
                    else:
                        result.append(currentPoint)

        return result

    def findBiPoints(self, line, y, x1, x2, tolerance):
        p1 = None

        for point in line:
            if abs(point.y - y) < tolerance:
                # x1 < x2, and we move from left to right
                if abs(point.x - x1) < tolerance:
                    p1 = point
                elif abs(point.x - x2) < tolerance:
                    return (p1, point)

        return (None, None)

    def perturbeAllData(self, lowData, medData, highData, lowZoneData, medZoneData, GUIMessage):
        lowData = self.perturbeData(lowData, self.lowROI.perturbeData(), GUIMessage + " perturbe lowData")
        medData = self.perturbeData(medData, self.medROI.perturbeData(), GUIMessage + " perturbe medData")
        highData = self.perturbeData(highData, self.highROI.perturbeData(), GUIMessage + " perturbe highData")

        lowZoneData = self.perturbeData(lowZoneData, self.lowROI.perturbeData(), GUIMessage + " perturbe lowZoneData")
        medZoneData = self.perturbeData(medZoneData, self.medROI.perturbeData(), GUIMessage + " perturbe medZoneData")

        return (lowData, medData, highData, lowZoneData, medZoneData)

