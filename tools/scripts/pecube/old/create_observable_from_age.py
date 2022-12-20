#!/usr/bin/env python

# written by Willi Kappler
# willi.kappler@uni-tuebingen.de

import os
import sys

print "generate observable file from ages file"

if len(sys.argv) != 2:
    print "you must provide a age file name"
    sys.exit(1)


input_file_name = sys.argv[1]

number_of_items = 0

all_ages = []

with open(input_file_name, "r") as input_file:
    for n, line in enumerate(input_file.readlines(), 1):
        if n == 4:
            key_value_pair = line.split(",")
            key, value = key_value_pair[0].split("=")
            if key.strip() == "n":
                number_of_items = int(value.strip())
                print "number_of_items: {0:d}".format(number_of_items)
        elif n > 4:
            all_ages.append(line.split())
            if n >= number_of_items + 4:
                break

print "number of items: {0:d}".format(len(all_ages))
print "first item:"
print all_ages[0]

with open("observable.dat", "w") as output_file:
    output_file.write("{0:d}\n".format(number_of_items))
    for values in all_ages:
        # x, y, z
        output_file.write("{0:s} {1:s} {2:s} 1 ".format(values[0], values[1], values[2]))
        output_file.write("{0:s} 0.5 ".format(values[3])) # AHe
        output_file.write("{0:s} 0.5 ".format(values[10])) # AFT
        output_file.write("{0:s} 0.5 ".format(values[12])) # ZHe
        output_file.write("{0:s} 0.5 ".format(values[13])) # ZFT
        output_file.write("{0:s} 0.5 ".format(values[17])) # MAr
        for entry in values[13:]:
            output_file.write("{0:s} 0.5 ".format(entry))
        output_file.write("\n")


