#!/usr/bin/env python

import sys
import os

print sys.argv

if len(sys.argv) == 1:
    print "no input files given!"
    sys.exit(1)

global_erate_counter = {}

with open("exported_good_fits_all.txt", "w") as output_file:
    output_file.write("# number_of_good_fits(j), k, j, chi_squared(k, j), predicted_ages(k, j)\n")
    for file_name in sys.argv[1:]:
        with open(file_name, "r") as input_file:
            print "opening file: {0:s}".format(file_name)
            local_erate_counter = {}
            for line in input_file.readlines():
                line_splitted = line.split()

                if len(line_splitted) == 5:
                    num_of_good_fits = int(line_splitted[0])
                    index_k = int(line_splitted[1])
                    index_j = int(line_splitted[2])
                    chi_squared = float(line_splitted[3])
                    predicted_age = float(line_splitted[4])

                    if (index_k, index_j) in local_erate_counter:
                        if num_of_good_fits > local_erate_counter[(index_k, index_j)]:
                            local_erate_counter[(index_k, index_j)] = num_of_good_fits
                    else:
                        local_erate_counter[(index_k, index_j)] = num_of_good_fits

                    if (index_k, index_j) in global_erate_counter:
                        output_file.write("{0:d} {1:d} {2:d} {3:f} {4:f}\n".format(num_of_good_fits +
                            global_erate_counter[(index_k, index_j)], index_k, index_j, chi_squared, predicted_age))
                    else:
                        output_file.write("{0:d} {1:d} {2:d} {3:f} {4:f}\n".format(num_of_good_fits,
                            index_k, index_j, chi_squared, predicted_age))

        for (index_k, index_j) in local_erate_counter.keys():
            if (index_k, index_j) in global_erate_counter:
                global_erate_counter[(index_k, index_j)] = global_erate_counter[(index_k, index_j)] + local_erate_counter[(index_k, index_j)]
            else:
                global_erate_counter[(index_k, index_j)] = local_erate_counter[(index_k, index_j)]

        print "local counter:"
        print local_erate_counter
        print "global counter:"
        print global_erate_counter
        print "\n----------\n"

