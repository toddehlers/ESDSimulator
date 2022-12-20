# written by Willi Kappler
# willi.kappler@uni-tuebingen.de
#

# system import
import os

# local import
import baseConfig
import boolrange
import stringrange
import intrange
import floatrange

class IceCascadeConfig(baseConfig.BaseConfig):
    def __init__(self, filename, version, numCPU):
        baseConfig.BaseConfig.__init__(self, filename, version, numCPU)

        self.outputDir = "icecascade"
        self.executable = "icecascade"

# The order of parameters is important!
# It has to match the order in the input file:

        self.parameters.append((boolrange.BoolRange("fluvial_erosion")))
        self.parameters.append((boolrange.BoolRange("ideposition")))
        self.parameters.append((boolrange.BoolRange("idiffusion")))
        self.parameters.append((boolrange.BoolRange("ilandslide")))
        self.parameters.append((boolrange.BoolRange("iceflag")))
        self.parameters.append((boolrange.BoolRange("glacial_erosion")))
        self.parameters.append((boolrange.BoolRange("iflexure")))
        self.parameters.append((boolrange.BoolRange("ihorizontal")))
        self.parameters.append((boolrange.BoolRange("iflag_oro")))
        self.parameters.append((boolrange.BoolRange("tecflag")))
        self.parameters.append((intrange.IntRange("nshortwrite")))
        self.parameters.append((intrange.IntRange("iflux")))
        self.parameters.append((intrange.IntRange("writetime")))
        self.parameters.append((intrange.IntRange("writeice")))
        self.parameters.append((stringrange.StringRange("run_name")))
        self.parameters.append((boolrange.BoolRange("iadjust")))
        self.parameters.append((floatrange.FloatRange("dt")))
        self.parameters.append((floatrange.FloatRange("endtime")))
        self.parameters.append((intrange.IntRange("nx")))
        self.parameters.append((intrange.IntRange("ny")))
        self.parameters.append((floatrange.FloatRange("sidex")))
        self.parameters.append((floatrange.FloatRange("sidey")))
        self.parameters.append((boolrange.BoolRange("iadapt")))
        self.parameters.append((boolrange.BoolRange("meshread")))
        self.parameters.append((stringrange.StringRange("meshname")))
        self.parameters.append((intrange.IntRange("meshformat")))
        self.parameters.append((boolrange.BoolRange("addshelf")))
        self.parameters.append((intrange.IntRange("Exr")))
        self.parameters.append((intrange.IntRange("Exl")))
        self.parameters.append((intrange.IntRange("Eyd")))
        self.parameters.append((intrange.IntRange("Eyu")))
        self.parameters.append((floatrange.FloatRange("slopexr")))
        self.parameters.append((floatrange.FloatRange("slopexl")))
        self.parameters.append((floatrange.FloatRange("slopeyd")))
        self.parameters.append((floatrange.FloatRange("slopyu")))
        self.parameters.append((intrange.IntRange("nxe")))
        self.parameters.append((intrange.IntRange("nye")))
        self.parameters.append((floatrange.FloatRange("xkf")))
        self.parameters.append((floatrange.FloatRange("xlf_BR")))
        self.parameters.append((floatrange.FloatRange("width_c")))
        self.parameters.append((floatrange.FloatRange("thresh")))
        self.parameters.append((floatrange.FloatRange("xlf_AL")))
        self.parameters.append((floatrange.FloatRange("sea_level")))
        self.parameters.append((floatrange.FloatRange("xkdiff")))
        self.parameters.append((intrange.IntRange("lsmeth")))
        self.parameters.append((floatrange.FloatRange("pmax")))
        self.parameters.append((floatrange.FloatRange("distmax")))
        self.parameters.append((floatrange.FloatRange("cohes")))
        self.parameters.append((floatrange.FloatRange("rho")))
        self.parameters.append((floatrange.FloatRange("grav")))
        self.parameters.append((floatrange.FloatRange("xk0")))
        self.parameters.append((floatrange.FloatRange("xk1")))
        self.parameters.append((floatrange.FloatRange("dtc")))
        self.parameters.append((floatrange.FloatRange("calc_rain_ice")))
        self.parameters.append((floatrange.FloatRange("shallow_ice")))
        self.parameters.append((floatrange.FloatRange("dh_allowed")))
        self.parameters.append((boolrange.BoolRange("itemp")))
        self.parameters.append((floatrange.FloatRange("ice_flow")))
        self.parameters.append((floatrange.FloatRange("sliding_law")))
        self.parameters.append((floatrange.FloatRange("power_law_exp")))
        self.parameters.append((floatrange.FloatRange("exp_sliding_law")))
        self.parameters.append((floatrange.FloatRange("rho_ice")))
        self.parameters.append((intrange.IntRange("constriction")))
        self.parameters.append((floatrange.FloatRange("erosion_rate")))
        self.parameters.append((floatrange.FloatRange("ice_ero_pow")))
        self.parameters.append((floatrange.FloatRange("basal_heat_flux")))
        self.parameters.append((floatrange.FloatRange("conductivity")))
        self.parameters.append((floatrange.FloatRange("critangle")))
        self.parameters.append((floatrange.FloatRange("calvecoef")))
        self.parameters.append((floatrange.FloatRange("uplift_rate")))
        self.parameters.append((floatrange.FloatRange("advec_velx")))
        self.parameters.append((floatrange.FloatRange("advec_vely")))
        self.parameters.append((floatrange.FloatRange("hflex")))
        self.parameters.append((boolrange.BoolRange("ixflex")))
        self.parameters.append((boolrange.BoolRange("iyflex")))
        self.parameters.append((floatrange.FloatRange("thickflex")))
        self.parameters.append((floatrange.FloatRange("young_mod")))
        self.parameters.append((floatrange.FloatRange("pratio")))
        self.parameters.append((floatrange.FloatRange("rhocflex")))
        self.parameters.append((floatrange.FloatRange("rhoaflex")))
        self.parameters.append((floatrange.FloatRange("flexon")))
        self.parameters.append((floatrange.FloatRange("calc_rain")))
        self.parameters.append((intrange.IntRange("del_gr")))
        self.parameters.append((floatrange.FloatRange("rain_vel")))
        self.parameters.append((intrange.IntRange("calc_ice")))
        self.parameters.append((floatrange.FloatRange("iceon")))
        self.parameters.append((floatrange.FloatRange("a0")))
        self.parameters.append((floatrange.FloatRange("a1")))
        self.parameters.append((floatrange.FloatRange("alf")))
        self.parameters.append((floatrange.FloatRange("wnd")))
        self.parameters.append((floatrange.FloatRange("angle")))
        self.parameters.append((floatrange.FloatRange("xwind_s")))
        self.parameters.append((floatrange.FloatRange("upwnd_s")))
        self.parameters.append((intrange.IntRange("subsam")))
        self.parameters.append((floatrange.FloatRange("xlapse_rate")))
        self.parameters.append((floatrange.FloatRange("At")))
        self.parameters.append((floatrange.FloatRange("xpmax")))
        self.parameters.append((boolrange.BoolRange("use_max_accu")))
        self.parameters.append((floatrange.FloatRange("kmelt")))
        self.parameters.append((floatrange.FloatRange("temp0min")))
        self.parameters.append((floatrange.FloatRange("temp0max")))
        self.parameters.append((floatrange.FloatRange("xp")))

    def loadConfig(self):
        inFile = open(self.filename, "r")

        for parameter in self.parameters:
            line = inFile.readline()
            line = line.strip()
            comments = []

            while line.startswith("#") or len(line) == 0:
                comments.append(line)
                line = inFile.readline()
                line = line.strip()

            parameter.comments = comments
            parameter.readLine(line)

        inFile.close()

    def writeConfigFile(self, outputFile):
        for parameter in self.parameters:
            for comment in parameter.comments:
                outputFile.write(comment + "\n")
            outputFile.write(parameter.valueAsString() + "\n")

    def writeParameters(self):
        self.stepCounter = self.stepCounter + 1
        print "stepCounter: ", self.stepCounter

        currentDirectory = "{0:s}/{1:04d}".format(self.outputDir, self.stepCounter)
        inputDirectory = "{0:s}/input/IceCascade/".format(currentDirectory)

        try:
            os.mkdir(currentDirectory)
            os.symlink("../../cascade_template/input", "{0:s}/input".format(currentDirectory))
            os.symlink("../../cascade_template/icecascade", "{0:s}/icecascade".format(currentDirectory))
            os.mkdir("{0:s}/output".format(currentDirectory))
            os.mkdir("{0:s}/output/IceCascade".format(currentDirectory))
        except OSError:
            pass # using existing directory

        self.currentInputFile = os.path.basename(self.filename)
        paramOutFile = open("{0:s}/{1:s}".format(inputDirectory, self.currentInputFile), "w")

        self.writeConfigFile(paramOutFile)

        paramOutFile.close()


