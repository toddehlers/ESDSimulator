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

class PecubeConfig(baseConfig.BaseConfig):
    def __init__(self, filename, version, numCPU):
        baseConfig.BaseConfig.__init__(self, filename, version, numCPU)

        self.outputDir = "pecube"
        self.executable = "pecube"

        # input 1
        self.parameters.append((stringrange.StringRange("run_name")))

        # input 2
        self.parameters.append((intrange.IntRange("num_of_topo_files")))

        # input 3
        self.parameters.append((intrange.IntRange("topo_file_flag")))

        # input 4
        self.parameters.append((stringrange.StringRange("topo_files", sizeVar="num_of_topo_files")))

        # input 5
        self.parameters.append((intrange.IntRange("coord_system")))

        # input 6
        self.parameters.append((intrange.IntRange("nx", numOfParamInLine=2)))
        self.parameters.append((intrange.IntRange("ny")))

        # input 7
        self.parameters.append((floatrange.FloatRange("longitude_spacing", numOfParamInLine=2)))
        self.parameters.append((floatrange.FloatRange("latitude spacing")))

        # input 8
        self.parameters.append((intrange.IntRange("nskip")))

        # input 9
        self.parameters.append((floatrange.FloatRange("geo_location_x", numOfParamInLine=2)))
        self.parameters.append((floatrange.FloatRange("geo_location_y")))

        # input 10
        self.parameters.append((intrange.IntRange("time_steps")))

        # input 11
        self.parameters.append((floatrange.FloatRange("erosional_time_scale")))

        # input 12
        self.parameters.append((floatrange.FloatRange("time_in_my", numOfParamInLine=13, sizeVar="time_steps", sizeVarMod=1)))
        self.parameters.append((floatrange.FloatRange("amplification_factor")))
        self.parameters.append((floatrange.FloatRange("vertical_offset_factor")))
        self.parameters.append((floatrange.FloatRange("time_temperature_flag")))
        self.parameters.append((floatrange.FloatRange("kinematic_flag")))
        self.parameters.append((floatrange.FloatRange("kinematic_details")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_g")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_h")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_i")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_j")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_k")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_l")))
        self.parameters.append((floatrange.FloatRange("kinematic_param_m")))

        # input 12b
        self.parameters.append((stringrange.StringRange("velocity_files", sizeVar="time_steps", sizeVarMod=1)))

        # input 13
        self.parameters.append((intrange.IntRange("isostasy_flag", numOfParamInLine=6)))
        self.parameters.append((floatrange.FloatRange("young_modules")))
        self.parameters.append((floatrange.FloatRange("poisson_ration")))
        self.parameters.append((floatrange.FloatRange("elastic_plate_thickness")))
        self.parameters.append((intrange.IntRange("fft_grid_x")))
        self.parameters.append((intrange.IntRange("fft_grid_y")))

        # input 14
        self.parameters.append((floatrange.FloatRange("model_thickness", numOfParamInLine=6)))
        self.parameters.append((intrange.IntRange("num_of_planes")))
        self.parameters.append((floatrange.FloatRange("thermal_condictivity")))
        self.parameters.append((floatrange.FloatRange("heat_capacity")))
        self.parameters.append((floatrange.FloatRange("crustal_density")))
        self.parameters.append((floatrange.FloatRange("mantle_density")))
        self.parameters.append((floatrange.FloatRange("base_temperature", numOfParamInLine=8)))
        self.parameters.append((floatrange.FloatRange("surface_temperature")))
        self.parameters.append((floatrange.FloatRange("atmo_lapse_rate")))
        self.parameters.append((floatrange.FloatRange("crustal_heat_production")))
        self.parameters.append((floatrange.FloatRange("e_folding_depth")))
        self.parameters.append((floatrange.FloatRange("mantle_heat_production")))
        self.parameters.append((intrange.IntRange("shear_heating_flag")))
        self.parameters.append((intrange.IntRange("shear_heating_constant")))

        # input 15
        self.parameters.append((floatrange.FloatRange("vol_heat_production_1", numOfParamInLine=5)))
        self.parameters.append((floatrange.FloatRange("vol_heat_production_2")))
        self.parameters.append((floatrange.FloatRange("vol_heat_production_3")))
        self.parameters.append((floatrange.FloatRange("vol_heat_production_4")))
        self.parameters.append((floatrange.FloatRange("vol_heat_production_5")))
        self.parameters.append((floatrange.FloatRange("thermal_conductivity_1", numOfParamInLine=5)))
        self.parameters.append((floatrange.FloatRange("thermal_conductivity_2")))
        self.parameters.append((floatrange.FloatRange("thermal_conductivity_3")))
        self.parameters.append((floatrange.FloatRange("thermal_conductivity_4")))
        self.parameters.append((floatrange.FloatRange("thermal_conductivity_5")))
        self.parameters.append((floatrange.FloatRange("rock_density_1", numOfParamInLine=5)))
        self.parameters.append((floatrange.FloatRange("rock_density_2")))
        self.parameters.append((floatrange.FloatRange("rock_density_3")))
        self.parameters.append((floatrange.FloatRange("rock_density_4")))
        self.parameters.append((floatrange.FloatRange("rock_density_5")))
        self.parameters.append((floatrange.FloatRange("spec_heat_capacity_1", numOfParamInLine=5)))
        self.parameters.append((floatrange.FloatRange("spec_heat_capacity_2")))
        self.parameters.append((floatrange.FloatRange("spec_heat_capacity_3")))
        self.parameters.append((floatrange.FloatRange("spec_heat_capacity_4")))
        self.parameters.append((floatrange.FloatRange("spec_heat_capacity_5")))

        # input 16
        self.parameters.append((intrange.IntRange("num_thermocron_data")))
        self.parameters.append((stringrange.StringRange("thermochron_data_file", sizeVar="num_thermocron_data")))

        # input 17
        self.parameters.append((intrange.IntRange("ages_output_flag_1", numOfParamInLine=11)))
        self.parameters.append((intrange.IntRange("ages_output_flag_2")))
        self.parameters.append((intrange.IntRange("ages_output_flag_3")))
        self.parameters.append((intrange.IntRange("ages_output_flag_4")))
        self.parameters.append((intrange.IntRange("ages_output_flag_5")))
        self.parameters.append((intrange.IntRange("ages_output_flag_6")))
        self.parameters.append((intrange.IntRange("ages_output_flag_7")))
        self.parameters.append((intrange.IntRange("ages_output_flag_8")))
        self.parameters.append((intrange.IntRange("ages_output_flag_9")))
        self.parameters.append((intrange.IntRange("ages_output_flag_10")))
        self.parameters.append((intrange.IntRange("ages_output_flag_11")))

        # input 18
        self.parameters.append((intrange.IntRange("detrital_age_flag")))

        # input 19
        self.parameters.append((intrange.IntRange("min_num_nodes")))

        # input 20
        self.parameters.append((stringrange.StringRange("cascade_tecplot_folder")))



    def loadConfig(self):
        inFile = open(self.filename, "r")
        paramCounter = 0

        def getLine(inFile, parameter):
            line = inFile.readline()
            line = line.strip()
            comments = []

            while line.startswith("$") or len(line) == 0:
                comments.append(line)
                line = inFile.readline()
                line = line.strip()

            parameter.comments = comments
            return line

        def expandParameter(parameter, paramCounter):
            sizeVarVal = None
            for p2 in self.parameters:
                if p2.name == parameter.sizeVar:
                    sizeVarVal = p2.startValue + parameter.sizeVarMod - 1
                    parameter.sizeVarVal = sizeVarVal
                    break

            if sizeVarVal == None:
                print "variable '{0:s}' not found for parameter '{1:s}'!".format(parameter.sizeVar, parameter.name)
                sys.exit(1)

            if sizeVarVal == -1:
#                print "sizeVarVal == -1 for '{0:s}'".format(parameter.name)
#                print "paramCounter: '{0:d}'".format(paramCounter)
#                print "name1: {0:s}".format(self.parameters[paramCounter].name)
                self.parameters.pop(paramCounter) # this parameter is not needed and can be removed
#                print "name2: {0:s}".format(self.parameters[paramCounter].name)

            newParams = []

            while sizeVarVal > 0:
                for i in xrange(parameter.numOfParamInLine):
                    newParams.append(self.parameters[paramCounter + i].copy())
                sizeVarVal = sizeVarVal - 1

#            print "sizeVar: '{0:s}', sizeVarVal: '{1:d}', numOfParamInLine: '{2:d}', new params: '{3:d}'".format(parameter.sizeVar, parameter.sizeVarVal, parameter.numOfParamInLine, len(newParams))

            insertIndex = paramCounter + parameter.numOfParamInLine

            self.parameters = self.parameters[:insertIndex] + newParams + self.parameters[insertIndex:]

        while paramCounter < len(self.parameters):
            parameter = self.parameters[paramCounter]

            if parameter.sizeVar:
                expandParameter(parameter, paramCounter)

            parameter = self.parameters[paramCounter]

            line = getLine(inFile, parameter)
            lineParam = line.split()

            actualNumOfParams = len(lineParam)

            if parameter.numOfParamInLine != actualNumOfParams:
                print "expected and actual number of parameters do not match for '{0:s}'\n".format(parameter.name)
                print "expected number of parameters in line: '{0:d}'".format(parameter.numOfParamInLine)
                print "actuall number of parameters in line: '{0:d}'".format(actualNumOfParams)
                print "current line: '{0:s}'".format(line)
                print parameter
                sys.exit(1)


            for param in lineParam:
                try:
                    self.parameters[paramCounter].readLine(param)
                except ValueError as error:
                    print "current line: '{0:s}'".format(line)
                    print "parameter: '{0:s}', paramCounter: '{1:d}', value: '{2:s}'".format(self.parameters[paramCounter].name, paramCounter, param)
                    raise error

                paramCounter = paramCounter + 1

        inFile.close()

    def writeConfigFile(self, outputFile):
        paramCounter = 0
        while paramCounter < len(self.parameters):
            parameter = self.parameters[paramCounter]
            for comment in parameter.comments:
                outputFile.write(comment + "\n")
            for i in xrange(parameter.numOfParamInLine):
                outputFile.write(self.parameters[paramCounter].valueAsString() + " ")
                paramCounter = paramCounter + 1
            outputFile.write("\n")

    def writeParameters(self):
        self.stepCounter = self.stepCounter + 1
        print "stepCounter: ", self.stepCounter

        currentDirectory = "{0:s}/{1:04d}".format(self.outputDir, self.stepCounter)

        try:
            os.mkdir(currentDirectory)
            os.symlink("../../pecube_template/pecube", "{0:s}/pecube".format(currentDirectory))
            os.symlink("../../pecube_template/input", "{0:s}/input".format(currentDirectory))
            os.mkdir("{0:s}/output".format(currentDirectory))
            os.mkdir("{0:s}/output/Pecube-D".format(currentDirectory))
        except OSError:
            pass # using existing directory

        self.currentInputFile = os.path.basename(self.filename)
        paramOutFile = open("{0:s}/{1:s}".format(currentDirectory, self.currentInputFile), "w")

        self.writeConfigFile(paramOutFile)

        paramOutFile.close()


