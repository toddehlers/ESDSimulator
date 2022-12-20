import Tkinter
import tkFileDialog
import tkMessageBox
import math

class CheckedEntry(Tkinter.Entry):
    def __init__(self, parent, textVar, width=None):
        Tkinter.Entry.__init__(self, parent, textvariable=textVar, width=width)
        self.textVar = textVar
        self.oldValue = None
        self.typeName = None
        self.bind("<FocusOut>", self.checkValue)
        self.bind("<FocusIn>", self.storeValue)

    def grid(self, row, column):
        Tkinter.Entry.grid(self, row=row, column=column, sticky="w", padx=2, pady=2)

    def storeValue(self, event):
        self.oldValue = self.textVar.get()

    def checkValue(self, event):
        valueAsString = self.textVar.get()

        try:
            valueAsFloat = self.convert(valueAsString)
        except ValueError:
            tkMessageBox.showerror("Invalid input", "This is not a valid\n{0:s} value: '{1:s}'".format(self.typeName, valueAsString))
            self.textVar.set(self.oldValue)

    def convert(self, value):
        raise NotImplementedError

class IntEntry(CheckedEntry):
    def __init__(self, parent, textVar, width):
        CheckedEntry.__init__(self, parent, textVar, width)
        self.oldValue = "0"
        self.typeName = "int"

    def convert(self, value):
        return int(value)

class FloatEntry(CheckedEntry):
    def __init__(self, parent, textVar, width=None):
        CheckedEntry.__init__(self, parent, textVar, width)
        self.oldValue = "0.0"
        self.typeName = "float"

    def convert(self, value):
        return float(value)

class ROIWidget(Tkinter.Frame):
    def __init__(self, parent, title, color, canvas, controller):
        Tkinter.Frame.__init__(self, parent, relief=Tkinter.GROOVE, borderwidth=2)
        self.title = title
        self.color = color
        self.canvas = canvas
        self.controller = controller
        self.smallerROI = None

        self.rubberBandID = None
        self.rubberBandBorder1ID = None
        self.rubberBandBorder2ID = None
        self.rubberBandX1 = None
        self.rubberBandY1 = None
        self.rubberBandX2 = None
        self.rubberBandY2 = None

        self.modelX = Tkinter.StringVar(value="0.0")
        self.modelY = Tkinter.StringVar(value="0.0")
        self.modelW = Tkinter.StringVar(value="0.0")
        self.modelH = Tkinter.StringVar(value="0.0")
        self.step = Tkinter.StringVar(value="10")
        self.perturbeDataCheck = Tkinter.IntVar(value=1)
        self.exportROICheck = Tkinter.IntVar(value=1)

        Tkinter.Label(self, text=self.title, bg=self.color, anchor="w").grid(row=0, column=0, columnspan=2, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="x:", anchor="e").grid(row=1, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="y:", anchor="e").grid(row=2, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="width:", anchor="e").grid(row=3, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="height:", anchor="e").grid(row=4, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="step:", anchor="e").grid(row=5, column=0, sticky="w", padx=2, pady=2)

        entry = FloatEntry(self, self.modelX, width=10)
        entry.grid(row=1, column=1)
        entry.bind("<FocusOut>", self.checkXW, add="+")
        entry.bind("<Return>", self.checkXW)

        entry = FloatEntry(self, self.modelY, width=10)
        entry.grid(row=2, column=1)
        entry.bind("<FocusOut>", self.checkYH, add="+")
        entry.bind("<Return>", self.checkYH)

        entry = FloatEntry(self, self.modelW, width=10)
        entry.grid(row=3, column=1)
        entry.bind("<FocusOut>", self.checkXW, add="+")
        entry.bind("<Return>", self.checkXW)

        entry = FloatEntry(self, self.modelH, width=10)
        entry.grid(row=4, column=1)
        entry.bind("<FocusOut>", self.checkYH, add="+")
        entry.bind("<Return>", self.checkYH)

        IntEntry(self, self.step, width=10).grid(row=5, column=1)

        Tkinter.Checkbutton(self, text="perturbe x,y pos.", variable=self.perturbeDataCheck, anchor="w").grid(row=6, column=0, columnspan=2, sticky="w", padx=2, pady=2)
        Tkinter.Checkbutton(self, text="export", variable=self.exportROICheck, anchor="w").grid(row=7, column=0, columnspan=2, sticky="w", padx=2, pady=2)

    def getStep(self):
        return int(self.step.get())

    def reset(self):
        self.step.set("10")
        self.clear()

    def checkXW(self, event):
        mx, mw = self.controller.checkXW(float(self.modelX.get()),float(self.modelW.get()))

        self.modelX.set(mx)
        self.modelW.set(mw)

        self.valueEntered()

    def checkYH(self, event):
        my, mh = self.controller.checkYH(float(self.modelY.get()),float(self.modelH.get()))

        self.modelY.set(my)
        self.modelH.set(mh)

        self.valueEntered()

    def valueEntered(self):
        if self.controller.modelLoaded:
            ROI = self.controller.convertModelToROI(*self.getModelROI())

            self.setRubberBandOrigin(ROI[0], ROI[1])
            self.updateRubberBand(ROI[2], ROI[3])

    def getModelROI(self):
        return (float(self.modelX.get()), float(self.modelY.get()), float(self.modelW.get()), float(self.modelH.get()))

    def setRubberBandOrigin(self, x1, y1):
        self.rubberBandX1 = x1
        self.rubberBandY1 = y1

    def updateRubberBand(self, x2, y2):
        if x2 < 0:
            x2 = 0
        if x2 > self.controller.screenWidth:
            x2 = self.controller.screenWidth

        if y2 < 0:
            y2 = 0
        if y2 > self.controller.screenHeight:
            y2 = self.controller.screenHeight

        self.rubberBandX2 = x2
        self.rubberBandY2 = y2

        self.canvas.delete(self.rubberBandID)
        self.canvas.delete(self.rubberBandBorder1ID)
        self.canvas.delete(self.rubberBandBorder2ID)
        self.rubberBandBorder1ID = self.canvas.create_rectangle(self.rubberBandX1-1, self.rubberBandY1-1, self.rubberBandX2+1, self.rubberBandY2+1, fill="", outline="#000000")
        self.rubberBandBorder2ID = self.canvas.create_rectangle(self.rubberBandX1+1, self.rubberBandY1+1, self.rubberBandX2-1, self.rubberBandY2-1, fill="", outline="#000000")
        self.rubberBandID = self.canvas.create_rectangle(self.rubberBandX1, self.rubberBandY1, self.rubberBandX2, self.rubberBandY2, fill="", outline=self.color)

    def updateTextValues(self):
        tempX1 = min(self.rubberBandX1, self.rubberBandX2)
        tempY1 = min(self.rubberBandY1, self.rubberBandY2)
        tempX2 = max(self.rubberBandX1, self.rubberBandX2)
        tempY2 = max(self.rubberBandY1, self.rubberBandY2)

        self.rubberBandX1 = tempX1
        self.rubberBandY1 = tempY1
        self.rubberBandX2 = tempX2
        self.rubberBandY2 = tempY2

        model = self.controller.convertROIToModel(self.rubberBandX1, self.rubberBandY1, self.rubberBandX2, self.rubberBandY2)

        self.modelX.set(model[0])
        self.modelY.set(model[1])
        self.modelW.set(model[2])
        self.modelH.set(model[3])

    def perturbeData(self):
        return self.perturbeDataCheck.get()

    def exportROI(self):
        return self.exportROICheck.get()

    def __lt__(self, other):
        if not self.exportROI():
            return False

        if other == None:
            return False

        if not other.exportROI():
            return self < other.smallerROI

        if self.rubberBandX1 > other.rubberBandX1:
            return True

        if self.rubberBandY1 > other.rubberBandY1:
            return True

        if self.rubberBandX2 < other.rubberBandX2:
            return True

        if self.rubberBandY2 < other.rubberBandY2:
            return True

        if self.getStep() < other.getStep():
            return True

        return False

    def getSettings(self):
        return (self.modelX.get(), self.modelY.get(), self.modelW.get(), self.modelH.get(), self.step.get(), self.perturbeDataCheck.get(), self.exportROICheck.get())

    def setSettings(self, settings):
        self.modelX.set(settings[0])
        self.modelY.set(settings[1])
        self.modelW.set(settings[2])
        self.modelH.set(settings[3])
        self.step.set(settings[4].strip())
        self.perturbeDataCheck.set(int(settings[5]))
        self.exportROICheck.set(int(settings[6]))

        self.checkXW(None)
        self.checkYH(None)

    def setCoordinatesToGrid(self, gridSpacing, outerGrid):
        x = float(self.modelX.get())
        y = float(self.modelY.get())
        w = float(self.modelW.get())
        h = float(self.modelH.get())

        print "setCoordinatesToGrid, old values: x: {0:f}, y: {1:f}, w: {2:f}, h: {3:f}".format(x, y, w, h)

        if outerGrid:
            self.modelX.set(math.floor(x/gridSpacing) * gridSpacing)
            self.modelY.set(math.floor(y/gridSpacing) * gridSpacing)
            self.modelW.set((math.ceil((x + w)/gridSpacing) * gridSpacing) - x)
            self.modelH.set((math.ceil((y + h)/gridSpacing) * gridSpacing) - y)
        else:
            self.modelX.set(math.ceil(x/gridSpacing) * gridSpacing)
            self.modelY.set(math.ceil(y/gridSpacing) * gridSpacing)
            self.modelW.set((math.floor((x + w)/gridSpacing) * gridSpacing) - x)
            self.modelH.set((math.floor((y + h)/gridSpacing) * gridSpacing) - y)

        print "setCoordinatesToGrid, new values: x: {0:s}, y: {1:s}, w: {2:s}, h: {3:s}\n".format(self.modelX.get(), self.modelY.get(), self.modelW.get(), self.modelH.get())

    def clear(self):
        self.modelX.set("0.0")
        self.modelY.set("0.0")
        self.modelW.set("0.0")
        self.modelH.set("0.0")
        self.canvas.delete(self.rubberBandID)
        self.canvas.delete(self.rubberBandBorder1ID)
        self.canvas.delete(self.rubberBandBorder2ID)
        self.rubberBandX1 = 0
        self.rubberBandY1 = 0
        self.rubberBandX2 = 0
        self.rubberBandY2 = 0

class ChildBorderWidget(Tkinter.Frame):
    def __init__(self, parent):
        Tkinter.Frame.__init__(self, parent, relief=Tkinter.GROOVE, borderwidth=2)

        self.borderTop = Tkinter.StringVar(value="open")
        self.borderBottom = Tkinter.StringVar(value="open")
        self.borderLeft = Tkinter.StringVar(value="open")
        self.borderRight = Tkinter.StringVar(value="open")
        self.borderValue = Tkinter.StringVar(value="0.0")
        self.fixedBorderCheck = Tkinter.IntVar(value=0)

        Tkinter.Label(self, text="Child border settings:", anchor="e").grid(row=0, column=0, columnspan=2, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="top:", anchor="e").grid(row=1, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="bottom:", anchor="e").grid(row=2, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="left:", anchor="e").grid(row=3, column=0, sticky="w", padx=2, pady=2)
        Tkinter.Label(self, text="right:", anchor="e").grid(row=4, column=0, sticky="w", padx=2, pady=2)

        # values for border: close = 1, open = 2
        Tkinter.OptionMenu(self, self.borderTop, "open", "close").grid(row=1, column=1, sticky="w", padx=2, pady=2)
        Tkinter.OptionMenu(self, self.borderBottom, "open", "close").grid(row=2, column=1, sticky="w", padx=2, pady=2)
        Tkinter.OptionMenu(self, self.borderLeft, "open", "close").grid(row=3, column=1, sticky="w", padx=2, pady=2)
        Tkinter.OptionMenu(self, self.borderRight, "open", "close").grid(row=4, column=1, sticky="w", padx=2, pady=2)

        Tkinter.Checkbutton(self, text="fixed value:", variable=self.fixedBorderCheck, anchor="w").grid(row=5, column=0, sticky="w", padx=2, pady=2)
        FloatEntry(self, self.borderValue, width=8).grid(row=5, column=1)

    def getValues(self):
        topValue = 2
        bottomValue = 2
        leftValue = 2
        rightValue = 2

        if self.borderTop.get() == "close":
            topValue = 1

        if self.borderBottom.get() == "close":
            bottomValue = 1

        if self.borderLeft.get() == "close":
            leftValue = 1

        if self.borderRight.get() == "close":
            rightValue = 1

        return (topValue, bottomValue, leftValue, rightValue)

    def setValues(self, topValue, bottomValue, leftValue, rightValue):
        def convertValue(value):
            if value == 1:
                return "close"
            else:
                return "open"

        self.borderTop.set(convertValue(topValue))
        self.borderBottom.set(convertValue(bottomValue))
        self.borderLeft.set(convertValue(leftValue))
        self.borderRight.set(convertValue(rightValue))

    def fixedBorder(self):
        if self.fixedBorderCheck.get():
            return float(self.borderValue.get())
        else:
            return None

class View(object):
    def __init__(self, controller):
        self.controller = controller

        self.currentROI = None

        # some constants
        self.LOW_MODE = 1
        self.LOW_MODE_COLOR = "#ff0000"
        self.MED_MODE = 2
        self.MED_MODE_COLOR = "#00ff00"
        self.HIGH_MODE = 3
        self.HIGH_MODE_COLOR = "#ffff00"
        self.ARBITRARY_MODE_COLOR = "#0000ff"

        self.startX = 0
        self.startY = 0
        self.endX = 0
        self.endY = 0

        self.root = Tkinter.Tk()
        self.root.title("convertArc")
        self.root.minsize(self.controller.screenWidth,self.controller.screenHeight)
        self.mainWindow = Tkinter.Frame(self.root, relief=Tkinter.GROOVE, borderwidth=4)
        self.mainWindow.pack(fill=Tkinter.BOTH, expand=Tkinter.YES)
        self.topFrame = Tkinter.Frame(self.mainWindow)
        self.topFrame.pack(side=Tkinter.TOP, fill=Tkinter.X, expand=Tkinter.NO)
        self.bottomFrame = Tkinter.Frame(self.mainWindow)
        self.bottomFrame.pack(side=Tkinter.TOP, fill=Tkinter.X, expand=Tkinter.NO)
        self.createWidgets()

        self.loadProgress = Tkinter.StringVar(value="Loading")

    def run(self):
        self.mainWindow.mainloop()
        self.root.destroy()

    def createWidgets(self):
        self.createMenu()
        self.createCanvas()
        self.createTopPanel()
        self.createRightPanel()

    def createMenu(self):
        menubar = Tkinter.Menu(self.root)
        fileMenu = Tkinter.Menu(menubar, tearoff=1)
        fileMenu.add_command(label="Open DEM...", command=self.fileMenuOpenDEM)
        fileMenu.add_command(label="Open GDAL file...", command=self.fileMenuOpenGDAL)
        fileMenu.add_command(label="Open XML...", command=self.fileMenuOpenXML)
        fileMenu.add_separator()
        fileMenu.add_command(label="Export Tecplot...", command=self.fileMenuExportTecplot)
        fileMenu.add_command(label="Export Cascade...", command=self.fileMenuExportCascade)
        fileMenu.add_command(label="Export Pecube...", command=self.fileMenuExportPecube)
        fileMenu.add_command(label="Export CSV...", command=self.fileMenuExportCSV)
        fileMenu.add_command(label="Export Child...", command=self.fileMenuExportChild)
        fileMenu.add_command(label="Export STL...", command=self.fileMenuExportSTL)
        fileMenu.add_command(label="Export ENVI Header...", command=self.fileMenuExportEnvi)
        fileMenu.add_separator()
        fileMenu.add_command(label="Save settings...", command=self.fileMenuSaveSettings)
        fileMenu.add_command(label="Load settings...", command=self.fileMenuLoadSettings)
        fileMenu.add_separator()
        fileMenu.add_command(label="Quit", command=self.fileMenuQuit)
        menubar.add_cascade(label="File", menu=fileMenu)

        editMenu = Tkinter.Menu(menubar, tearoff=1)
        editMenu.add_command(label="Set DEM origin", command=self.editMenuDEMOrigin)
        editMenu.add_command(label="Set export origin", command=self.editMenuExportOrigin)
        editMenu.add_separator()
        editMenu.add_command(label="Select low res ROI", command=self.editMenuLowROI)
        editMenu.add_command(label="Select med res ROI", command=self.editMenuMedROI)
        editMenu.add_command(label="Select high res ROI", command=self.editMenuHighROI)
        editMenu.add_command(label="Select arbitrary box", command=self.editMenuArbitraryBox)
        editMenu.add_command(label="Clear all ROI", command=self.editMenuClearAllROI)
        menubar.add_cascade(label="Edit", menu=editMenu)

        self.root.config(menu=menubar)

    def createCanvas(self):
        self.pic = Tkinter.PhotoImage(name="data", width=self.controller.screenWidth, height=self.controller.screenHeight)
        self.clearImage()

        self.canvas = Tkinter.Canvas(self.bottomFrame, width=self.controller.screenWidth, height=self.controller.screenHeight, relief=Tkinter.GROOVE, borderwidth=2)
        self.canvas.create_image((self.controller.screenWidth / 2,self.controller.screenHeight / 2), image=self.pic)
        self.canvas.pack(side=Tkinter.LEFT, fill=Tkinter.Y, expand=Tkinter.NO)
        self.canvas.bind("<Button-1>", self.leftMBClick)
        self.canvas.bind("<B1-Motion>", self.leftMBMotion)
        self.canvas.bind("<ButtonRelease-1>", self.leftMBRelease)

    def clearImage(self):
        imageList = []
        for i in xrange(self.controller.screenWidth):
            imageList.append("{")
            for j in xrange(self.controller.screenHeight):
                imageList.append("#000000")
            imageList.append("}")
        imageString = " ".join(imageList)
        self.pic.put(imageString)
        del imageList
        del imageString

    def createTopPanel(self):
        self.DEMFile = Tkinter.StringVar(value="no file loaded")
        self.XMLFile = Tkinter.StringVar(value="no file loaded")
        self.DEMOrigin = Tkinter.StringVar(value="(0, 0)")
        self.exportOrigin = Tkinter.StringVar(value="(0, 0)")
        self.DEMRows = Tkinter.StringVar(value="0")
        self.DEMCols = Tkinter.StringVar(value="0")
        self.DEMCellSize = Tkinter.StringVar(value="0")
        self.DEMArea = Tkinter.StringVar(value="(0, 0)")

        # row 0
        Tkinter.Label(self.topFrame, text="DEM origin:", anchor="e").grid(row=0, column=0, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMOrigin).grid(row=0, column=1, sticky="w", padx=2, pady=2)

        Tkinter.Frame(self.topFrame, width=30).grid(row=0, column=2) # separator

        Tkinter.Label(self.topFrame, text="DEM rows:", anchor="e").grid(row=0, column=3, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMRows).grid(row=0, column=4, sticky="w", padx=2, pady=2)

        Tkinter.Frame(self.topFrame, width=30).grid(row=0, column=5) # separator

        Tkinter.Label(self.topFrame, text="DEM Cell size:", anchor="e").grid(row=0, column=6, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMCellSize).grid(row=0, column=7, sticky="w", padx=2, pady=2)

        Tkinter.Frame(self.topFrame, width=30).grid(row=0, column=8) # separator

        Tkinter.Label(self.topFrame, text="DEM file:", anchor="e").grid(row=0, column=9, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMFile).grid(row=0, column=10, sticky="w", padx=2, pady=2)

        # row 1
        Tkinter.Label(self.topFrame, text="Export origin:", anchor="e").grid(row=1, column=0, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.exportOrigin).grid(row=1, column=1, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, text="DEM cols:", anchor="e").grid(row=1, column=3, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMCols).grid(row=1, column=4, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, text="DEM area (km):", anchor="e").grid(row=1, column=6, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.DEMArea).grid(row=1, column=7, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, text="XML metadata file:", anchor="e").grid(row=1, column=9, sticky="w", padx=2, pady=2)

        Tkinter.Label(self.topFrame, anchor="w", bg="white", textvariable=self.XMLFile).grid(row=1, column=10, sticky="w", padx=2, pady=2)

    def createRightPanel(self):
        self.zScaleFactor = Tkinter.StringVar(value="1.0")
        self.noDataValue = Tkinter.StringVar(value="0.0")
        self.DEMInUTM = Tkinter.IntVar(value=1)
        self.perturbeFactor = Tkinter.StringVar(value="0.5")
        self.geoTIFFBand = Tkinter.IntVar(value=1)

        leftFrame = Tkinter.Frame(self.bottomFrame)
        leftFrame.pack(side=Tkinter.LEFT, fill=Tkinter.Y, expand=Tkinter.NO)

        self.lowROI = ROIWidget(leftFrame, "low resolution ROI (km)", self.LOW_MODE_COLOR, self.canvas, self.controller)
        self.lowROI.grid(row=0, column=0, padx=5, pady=5, rowspan=4, sticky="w")

        self.medROI = ROIWidget(leftFrame, "med resolution ROI (km)", self.MED_MODE_COLOR, self.canvas, self.controller)
        self.medROI.grid(row=4, column=0, padx=5, pady=5, rowspan=4, sticky="w")

        self.highROI = ROIWidget(leftFrame, "high resolution ROI (km)", self.HIGH_MODE_COLOR, self.canvas, self.controller)
        self.highROI.grid(row=8, column=0, padx=5, pady=5, rowspan=4, sticky="w")

        self.arbitraryROI = ROIWidget(leftFrame, "arbitrary ROI (km)", self.ARBITRARY_MODE_COLOR, self.canvas, self.controller)

        self.lowROI.smallerROI = self.medROI
        self.medROI.smallerROI = self.highROI
        self.highROI.smallerROI = None
        self.arbitraryROI.smallerROI = None

        self.controller.setROI(self.lowROI, self.medROI, self.highROI)

        Tkinter.Label(leftFrame, text="no data value:", anchor="e").grid(row=0, column=1, sticky="w", padx=2, pady=2)
        FloatEntry(leftFrame, self.noDataValue, width=5).grid(row=0, column=2)

        Tkinter.Label(leftFrame, text="elev. scale factor: z *", anchor="e").grid(row=1, column=1, sticky="w", padx=2, pady=2)
        FloatEntry(leftFrame, self.zScaleFactor, width=5).grid(row=1, column=2)

        Tkinter.Checkbutton(leftFrame, text="input DEM in UTM", variable=self.DEMInUTM, anchor="w").grid(row=2, column=1, sticky="w", padx=2, pady=2)

        Tkinter.Label(leftFrame, text="pertubation factor:", anchor="e").grid(row=3, column=1, sticky="w", padx=2, pady=2)
        FloatEntry(leftFrame, self.perturbeFactor, width=5).grid(row=3, column=2)

        Tkinter.Label(leftFrame, text="geoTIFF band:", anchor="e").grid(row=4, column=1, sticky="w", padx=2, pady=2)
        IntEntry(leftFrame, self.geoTIFFBand, width=5).grid(row=4, column=2)

        self.childBorder = ChildBorderWidget(leftFrame)
        self.childBorder.grid(row=5, column=1, padx=5, pady=5, columnspan=2, sticky="w")

    def resetValues(self):
        self.DEMFile.set("")
        self.XMLFile.set("")
        self.DEMOrigin.set("(0,0)")
        self.exportOrigin.set("(0,0)")
        self.DEMRows.set("")
        self.DEMCols.set("")
        self.DEMCellSize.set("")
        self.DEMArea.set("")
        self.lowROI.reset()
        self.medROI.reset()
        self.highROI.reset()

    def fileMenuOpenDEM(self):
        selectedFile = tkFileDialog.askopenfilename(parent=self.root, title='Please select a DEM file')
        if len(selectedFile) > 0:
            try:
                self.clearImage()
                self.resetValues()
                self.DEMFile.set(selectedFile)
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.loadDEM(selectedFile)
                progressWindow.destroy()
            except ValueError as ve:
                tkMessageBox.showerror("Error reading DEM file", "An error occurred while reading the DEM file:\n{0:s}".format(str(ve)))
                self.controller.modelLoaded = False

    def fileMenuOpenGDAL(self):
        selectedFile = tkFileDialog.askopenfilename(parent=self.root, title='Please select a GDAL file')
        if len(selectedFile) > 0:
            try:
                self.clearImage()
                self.resetValues()
                self.DEMFile.set(selectedFile)
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.loadGDAL(selectedFile)
                progressWindow.destroy()
            except ValueError as ve:
                tkMessageBox.showerror("Error reading GDAL file", "An error occurred while reading the GDAL file:\n{0:s}".format(str(ve)))
                self.controller.modelLoaded = False

    def fileMenuOpenXML(self):
        self.notImplementedYet()
#        selectedFile = tkFileDialog.askopenfilename(parent=self.root, title='Please select a XML file')
#        if selectedFile != None:
#            self.XMLFile.set(selectedFile)
#            self.controller.loadXML(selectedFile)

    def fileMenuExportTecplot(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the tecplot output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportTecplot(selectedFile)
                progressWindow.destroy()

    def fileMenuExportCascade(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the cascade output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportCascade(selectedFile)
                progressWindow.destroy()

    def fileMenuExportPecube(self):
        selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the pecube output file')
        if len(selectedFile) > 0:
            progressWindow = Tkinter.Toplevel()
            Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
            self.controller.exportPecube(selectedFile)
            progressWindow.destroy()

    def fileMenuExportCSV(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the csv output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportCSV(selectedFile)
                progressWindow.destroy()

    def fileMenuExportChild(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the child output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportChild(selectedFile, self.childBorder.getValues(), self.childBorder.fixedBorder())
                progressWindow.destroy()

    def fileMenuExportSTL(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the STL output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportSTL(selectedFile, self.childBorder.getValues(), self.childBorder.fixedBorder())
                progressWindow.destroy()

    def fileMenuExportEnvi(self):
        if self.checkROI():
            selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the ENVI output file')
            if len(selectedFile) > 0:
                progressWindow = Tkinter.Toplevel()
                Tkinter.Label(progressWindow, textvariable=self.loadProgress).pack(ipadx=20, ipady=20)
                self.controller.exportEnvi(selectedFile, self.childBorder.getValues(), self.childBorder.fixedBorder())
                progressWindow.destroy()
    def fileMenuSaveSettings(self):
        selectedFile = tkFileDialog.asksaveasfilename(parent=self.root, title='Please select the settings file')
        if len(selectedFile) > 0:
            settingsFile = open(selectedFile, "w")
            settingsFile.write("# ROI settings: x, y, w, h, step, perturbe, export\n")
            settingsFile.write("LowROI: {0:s}, {1:s}, {2:s}, {3:s}, {4:s}, {5:d}, {6:d}\n".format(*self.lowROI.getSettings()))
            settingsFile.write("MedROI: {0:s}, {1:s}, {2:s}, {3:s}, {4:s}, {5:d}, {6:d}\n".format(*self.medROI.getSettings()))
            settingsFile.write("HighROI: {0:s}, {1:s}, {2:s}, {3:s}, {4:s}, {5:d}, {6:d}\n".format(*self.highROI.getSettings()))
            settingsFile.write("\n# Global settings: no data value, elevation scale factor, input in UTM, perturbation factor, geoTIFF band\n")
            settingsFile.write("Global: {0:s}, {1:s}, {2:d}, {3:s}, {4:d}\n".format(self.noDataValue.get(), self.zScaleFactor.get(), self.DEMInUTM.get(), self.perturbeFactor.get(), self.geoTIFFBand.get()))
            settingsFile.write("\n# Child border settings:  top, bottom, left, right\n")
            settingsFile.write("Border: {0:d}, {1:d}, {2:d}, {3:d}\n".format(*self.childBorder.getValues()))
            settingsFile.close()

    def fileMenuLoadSettings(self):
        selectedFile = tkFileDialog.askopenfilename(parent=self.root, title='Please select the settings file')
        if len(selectedFile) > 0:
            settingsFile = open(selectedFile)

            while True:
                line = settingsFile.readline()

                if len(line) == 0:
                    break

                if line.startswith("#"):
                    continue

                line = line.split(":")

                if len(line) == 1:
                    continue

                item = line[0].strip()
                values = line[1].split(",")

                if len(values) == 1:
                    continue

                if item == "LowROI":
                    self.lowROI.setSettings(values)
                elif item == "MedROI":
                    self.medROI.setSettings(values)
                elif item == "HighROI":
                    self.highROI.setSettings(values)
                elif item == "Global":
                    self.noDataValue.set(values[0].strip())
                    self.zScaleFactor.set(values[1].strip())
                    self.DEMInUTM.set(values[2].strip())
                    self.perturbeFactor.set(values[3].strip())
                    self.geoTIFFBand.set(values[4].strip())
                elif item == "Border":
                    self.childBorder.setValues(int(values[0]), int(values[1]), int(values[2]), int(values[3]))

            settingsFile.close()


    def fileMenuQuit(self):
        self.controller.quit()
        self.root.quit()

    def editMenuDEMOrigin(self):
        self.notImplementedYet()

    def editMenuExportOrigin(self):
        self.notImplementedYet()

    def editMenuLowROI(self):
        if self.controller.modelLoaded:
            self.currentROI = self.lowROI
            self.root.configure(cursor="crosshair")
        else:
            tkMessageBox.showerror("No DEM loaded", "You have to load a DEM file first!")


    def editMenuMedROI(self):
        if self.controller.modelLoaded:
            self.currentROI = self.medROI
            self.root.configure(cursor="crosshair")
        else:
            tkMessageBox.showerror("No DEM loaded", "You have to load a DEM file first!")

    def editMenuHighROI(self):
        if self.controller.modelLoaded:
            self.currentROI = self.highROI
            self.root.configure(cursor="crosshair")
        else:
            tkMessageBox.showerror("No DEM loaded", "You have to load a DEM file first!")

    def editMenuArbitraryBox(self):
        if self.controller.modelLoaded:
            self.currentROI = self.arbitraryROI
            self.root.configure(cursor="crosshair")
        else:
            tkMessageBox.showerror("No DEM loaded", "You have to load a DEM file first!")

    def editMenuClearAllROI(self):
        self.lowROI.clear()
        self.medROI.clear()
        self.highROI.clear()

    def leftMBClick(self, event):
        if self.currentROI != None:
            self.currentROI.setRubberBandOrigin(event.x, event.y)

    def leftMBMotion(self, event):
        if self.currentROI != None:
            self.currentROI.updateRubberBand(event.x, event.y)

    def leftMBRelease(self, event):
        if self.currentROI != None:
            self.currentROI.updateTextValues()
            self.currentROI = None
            self.root.configure(cursor="left_ptr")

    def setDEMProperties(self, ncols, nrows, cellsize, width, height):
        self.DEMCols.set(ncols)
        self.DEMRows.set(nrows)
        self.DEMCellSize.set(cellsize)
        self.DEMArea.set("({0:f}, {1:f})".format(width, height))

    def isUTM(self):
        return self.DEMInUTM.get()

    def getNODATA_replacement(self):
        return float(self.noDataValue.get())

    def getZFactor(self):
        return float(self.zScaleFactor.get())

    def getPertubeFactor(self):
        return float(self.perturbeFactor.get())

    def getGeoTIFFBand(self):
        return self.geoTIFFBand.get()

    def setDisplayData(self, data, minValue, maxValue):
        valueRange = maxValue - minValue
        result = []

        blueCyan = 0.25
        cyanGreen = 0.5
        greenYellow = 0.75
        colorFactor = 4.0 * 255.0

        for line in data:
            result.append("{")
            for value in line:
                value = (value - minValue) / valueRange
                if value <= blueCyan:
                    result.append("#00{0:02x}ff".format(int(value * colorFactor)))
                elif value > blueCyan and value <= cyanGreen:
                    result.append("#00ff{0:02x}".format(255 - int((value - blueCyan) * colorFactor)))
                elif value > cyanGreen and value <= greenYellow:
                    result.append("#{0:02x}ff00".format(int((value - cyanGreen) * colorFactor)))
                elif value > greenYellow:
                    result.append("#ff{0:02x}00".format(255 - int((value - greenYellow) * colorFactor)))

            result.append("}")

        self.pic.put(" ".join(result))

    def updateDisplay(self, currentValue, endValue, message=""):
        self.loadProgress.set("Processing ({0:d} of {1:d}) {2:s}".format(currentValue, endValue, message))
        self.mainWindow.update()

    def notImplementedYet(self):
        tkMessageBox.showerror("Not implemented yet", "This feature is not implemeted yet!")

    def checkROI(self):
        if self.lowROI < self.medROI:
            tkMessageBox.showerror("ROI error", "low resolution ROI must be outside of medium and outside of high resolution ROI!\nLow step must be bigger than medium or high step!")
            return False
        if self.medROI < self.highROI:
            tkMessageBox.showerror("ROI error", "medium resolution ROI must be outside of high resolution ROI!\nMedium step must be bigger than high step!")
            return False
        return True

