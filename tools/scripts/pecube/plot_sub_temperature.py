#!/usr/bin/env python3


import os.path
import sys
import struct

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

def process_file(temperature_file_name, number_of_nodes):
    print("temperature: {}".format(temperature_file_name))

    pos_x = []
    pos_z = []
    temp_values = []

    with open(temperature_file_name, "rb") as input_file:
        (num_of_sub_steps, current_step) = struct.unpack("=ii", input_file.read(4 + 4))
        print("current step: {}, number of sub steps: {}".format(current_step, num_of_sub_steps))

        for sub_step1 in range(1, num_of_sub_steps + 1):
            (dt, sub_step2, time_value) = struct.unpack("=did", input_file.read(8 + 4 + 8))
            if (sub_step1 != sub_step2):
                print("sub steps do not match:")
                print("sub step 1: {}, sub step 2: {}".format(sub_step1, sub_step2))
                sys.exit(1)
            for i in range(1, number_of_nodes + 1):
                (n_id, px, py, pz, temp) = struct.unpack("=idddd", input_file.read(4 + 8 + 8 + 8 + 8))
                if (i != n_id):
                    print("node id does not match:")
                    print("i: {}, n_id: {}".format(i, n_id))
                    sys.exit(1)
                if (py == 0.0) and (sub_step1 == 1):
                    pos_x.append(px)
                    pos_z.append(pz)
                    temp_values.append(temp)


    #time_steps.reverse()

    fig, all_axes = plt.subplots(1, 1, figsize=(20,10))
    all_axes.set_title("file: {}".format(temperature_file_name))

    sc = all_axes.scatter(pos_x, pos_z, c=temp_values, edgecolors="none")
    all_axes.grid(True)
    all_axes.set_xlabel("x [km]")
    all_axes.set_ylabel("z [km]")
    plt.colorbar(sc)

    figure_file = temperature_file_name.replace(".dat", ".png")

    fig.savefig(figure_file, bbox_inches="tight", dpi=300)
    plt.close(fig)

if __name__ == "__main__":

    number_of_nodes = int(sys.argv[1])

    file_index = 1

    while True:
        temperature_file_name = "temperature_field_sub_{0:04d}.dat".format(file_index)

        if os.path.exists(temperature_file_name):
            process_file(temperature_file_name, number_of_nodes)
            file_index += 1
        else:
            print("Finished!")
            break
