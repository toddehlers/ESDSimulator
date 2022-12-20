#!/bin/env python3


import os.path
import sys
import struct
import operator

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

def process_file(surface_file_name, borehole_file_name, file_index, point_of_interest):
    print("\n\nprocessing: '{}' and '{}'".format(surface_file_name, borehole_file_name))

    time_history_values1 = []
    time_history_values2 = []
    temperature_history_values1 = []
    temperature_history_values2 = []

    with open(surface_file_name, "rb") as input_file1:
        (num_of_points1, ) = struct.unpack("=i", input_file1.read(4))
        #print("num_of_points1: {}".format(num_of_points1))
        while True:
            (current_step1, ntime1) = struct.unpack("=ii", input_file1.read(4 + 4))
            #print("  current_step1: {}, ntime1: {}".format(current_step1, ntime1))
            for sub_step in range(1, ntime1 + 1):
                (sub_step1, time_history_value1) = struct.unpack("=id", input_file1.read(4 + 8))
                time_history_values1.append(time_history_value1)

            for sub_step in range(1, ntime1 + 1):
                #print("    sub_step: {}".format(sub_step))
                for current_id in range(1, num_of_points1 + 1):
                    #print("      current_id: {}".format(current_id))
                    (id1, temperature1) = struct.unpack("=id", input_file1.read(4 + 8))
                    #print("      id1: {}, temperature1: {}".format(id1, temperature1))
                    struct.unpack("=dddddd", input_file1.read(8 + 8 + 8 + 8 + 8 + 8))
                    if current_id == point_of_interest:
                        temperature_history_values1.append(temperature1)

            if current_step1 == 1:
                break

    combined = zip(time_history_values1, temperature_history_values1)
    combined = sorted(combined, key=operator.itemgetter(0))
    (time_history_values1, temperature_history_values1) = map(list, zip(*combined))

    print("read next file")

    with open(borehole_file_name, "rb") as input_file2:
        (num_of_points2, ) = struct.unpack("=i", input_file2.read(4))
        #print("num_of_points2: {}".format(num_of_points2))
        while True:
            (current_step2, ntime2) = struct.unpack("=ii", input_file2.read(4 + 4))
            #print("  current_step2: {}, ntime2: {}".format(current_step2, ntime2))
            for sub_step in range(1, ntime2 + 1):
                (sub_step2, time_history_value2) = struct.unpack("=id", input_file2.read(4 + 8))
                time_history_values2.append(time_history_value2)

            for sub_step in range(1, ntime2 + 1):
                #print("    sub_step: {}".format(sub_step))
                for current_id in range(1, num_of_points2 + 1):
                    #print("      current_id: {}".format(current_id))
                    (id2, temperature2) = struct.unpack("=id", input_file2.read(4 + 8))
                    #print("      id2: {}, temperature2: {}".format(id2, temperature2))
                    struct.unpack("=dddddd", input_file2.read(8 + 8 + 8 + 8 + 8 + 8))
                    if current_id == point_of_interest:
                        temperature_history_values2.append(temperature2)

            if current_step2 == 1:
                break

    combined = zip(time_history_values2, temperature_history_values2)
    combined = sorted(combined, key=operator.itemgetter(0))
    (time_history_values2, temperature_history_values2) = map(list, zip(*combined))

    output_file = "temperature_history_{}.png".format(file_index)
    fig, all_axes = plt.subplots(1, 1, figsize=(20,10))
    all_axes.set_title("file: {}".format(output_file))

    all_axes.plot(time_history_values1, temperature_history_values1, "b-", label="surface points")
    all_axes.plot(time_history_values2, temperature_history_values2, "g--", label="borehole points")
    all_axes.xaxis.set_major_locator(mti.MultipleLocator(1))
    all_axes.yaxis.set_major_locator(mti.MultipleLocator(10))
    all_axes.grid(True)
    all_axes.set_xlabel("Time [Ma]")
    all_axes.set_ylabel("Temperature [°C]")
    all_axes.set_xlim([0,22])
    all_axes.set_ylim([0,250])
    all_axes.legend(loc="upper right")
    all_axes.invert_xaxis()

    fig.savefig(output_file, bbox_inches="tight", dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    file_index = 14

    point_of_interest = int(sys.argv[1])

    while True:
        surface_file_name = "time_temperature_history_{0:04d}.bin".format(file_index)
        borehole_file_name = "borehole_time_temperature_history_{0:04d}.bin".format(file_index)

        if os.path.exists(surface_file_name) and os.path.exists(borehole_file_name):
            process_file(surface_file_name, borehole_file_name, file_index, point_of_interest)
            file_index += 1
        else:
            print("Finished!")
            break
