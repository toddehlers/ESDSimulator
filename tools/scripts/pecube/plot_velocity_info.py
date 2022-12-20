#!/bin/env python3

# 1, 65, 66, 67, 716, 2012, 2013, 2014, 2062, 2063

import os.path
import sys
import struct
import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

def process_file(velocity_file_name, temperature_file_name, surface_node):
    print("id: {}, velocity: {}, temperature: {}".format(surface_node, velocity_file_name, temperature_file_name))

    time_steps = []
    pos_x = []
    pos_z = []
    velo_x = []
    velo_z = []
    temperature = []
    time_values = []
    number_of_surface_nodes = 0
    start_step = 0

    with open(velocity_file_name, "rb") as input_file:
        (start_step, number_of_surface_nodes) = struct.unpack("ii", input_file.read(8))

        for current_step2 in range(start_step, 0, -1):
            (current_step1, ) = struct.unpack("i", input_file.read(4))
            if current_step1 != current_step2:
                    print("step does not match:")
                    print("current_step1: {} != current_step2: {}".format(current_step1, current_step2))
                    sys.exit(1)
            time_steps.append(current_step1)
            for current_node1 in range(1, number_of_surface_nodes + 1):
                (current_node2, ) = struct.unpack("i", input_file.read(4))
                if current_node1 != current_node2:
                    print("node id does not match:")
                    print("current_node1: {} != current_node2: {}".format(current_node1, current_node2))
                    sys.exit(1)
                (px, py, pz, vx, vy, vz) = struct.unpack("dddddd", input_file.read(8 * 6))
                if current_node1 == surface_node:
                    pos_x.append(px)
                    pos_z.append(pz)
                    velo_x.append(vx)
                    velo_z.append(vz)

    with open(temperature_file_name, "rb") as input_file:
        for current_step1 in range(start_step, 0, -1):
            (current_step2, number_of_sub_steps) = struct.unpack("ii", input_file.read(4 * 2))
            if current_step1 != current_step2:
                print("steps do not match:")
                print("step1: {}, step2: {}".format(current_step1, current_step2))
                sys.exit(1)
            #print("current_step: {}, number_of_sub_steps: {}".format(current_step1, number_of_sub_steps))
            for sub_step1 in range(1, number_of_sub_steps + 1):
                (sub_step2, time) = struct.unpack("=id", input_file.read(4 + 8))
                if sub_step1 != sub_step2:
                    print("sub steps do not match:")
                    print("sub_step1: {}, sub_step2: {}".format(sub_step1, sub_step2))
                    sys.exit(1)
                if sub_step1 == 1:
                    time_values.append(time)
                #print("sub_step: {}, time: {}".format(sub_step1, time))

            for sub_step1 in range(1, number_of_sub_steps + 1):
                for current_node1 in range(1, number_of_surface_nodes + 1):
                    (current_node2, ) = struct.unpack("i", input_file.read(4))
                    if current_node1 != current_node2:
                        print("node id does not match:")
                        print("current_node1: {} != current_node2: {}".format(current_node1, current_node2))
                        sys.exit(1)
                    (temp, ) = struct.unpack("d", input_file.read(8))
                    if (current_node1 == surface_node) and (sub_step1 == 1):
                        temperature.append(temp)


    #time_steps.reverse()

#    for (t1, t2, t3, px, pz) in zip(time_steps, time_values, temperature, pos_x, pos_z):
#        print("step: {}, time: {}, temperature: {}, px: {}, pz: {}".format(t1, t2, t3, px, pz))

    pos_x_min = min(pos_x)
    pos_x_max = max(pos_x)
    pos_z_min = min(pos_z)
    pos_z_max = max(pos_z)
    velo_x_min = min(velo_x)
    velo_x_max = max(velo_x)
    velo_z_min = min(velo_z)
    velo_z_max = max(velo_z)
    temperature_min = min(temperature)
    temperature_max = max(temperature)
    time_values_min = min(time_values)
    time_values_max = max(time_values)
    start_step_max = start_step

    if os.path.exists("velocity_info_scale_{}.txt".format(surface_node)):
        with open("velocity_info_scale_{}.txt".format(surface_node), "r") as f:
            values = [float(v) for v in f.readline().split()]
            pos_x_min = min(pos_x_min, values[0])
            pos_x_max = max(pos_x_max, values[1])
            pos_z_min = min(pos_z_min, values[2])
            pos_z_max = max(pos_z_max, values[3])
            velo_x_min = min(velo_x_min, values[4])
            velo_x_max = max(velo_x_max, values[5])
            velo_z_min = min(velo_z_min, values[6])
            velo_z_max = max(velo_z_max, values[7])
            temperature_min = min(temperature_min, values[8])
            temperature_max = max(temperature_max, values[9])
            time_values_min = min(time_values_min, values[10])
            time_values_max = max(time_values_max, values[11])
            start_step_max = max(start_step, values[12])

    with open("velocity_info_scale_{}.txt".format(surface_node), "w") as f:
        f.write("{} ".format(pos_x_min))
        f.write("{} ".format(pos_x_max))
        f.write("{} ".format(pos_z_min))
        f.write("{} ".format(pos_z_max))
        f.write("{} ".format(velo_x_min))
        f.write("{} ".format(velo_x_max))
        f.write("{} ".format(velo_z_min))
        f.write("{} ".format(velo_z_max))
        f.write("{} ".format(temperature_min))
        f.write("{} ".format(temperature_max))
        f.write("{} ".format(time_values_min))
        f.write("{} ".format(time_values_max))
        f.write("{} ".format(start_step_max))

    fig, all_axes = plt.subplots(6, 1, figsize=(10,20))
    all_axes[0].set_title("id: {}, file: {}".format(surface_node, velocity_file_name))

    all_axes[0].plot(time_steps, pos_x, '-bo')
    all_axes[0].grid(True)
    all_axes[0].set_xlabel("time step")
    all_axes[0].set_ylabel("position x [km]")
    all_axes[0].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[0].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[0].set_xlim(1, start_step_max)
    all_axes[0].set_ylim(pos_x_min, pos_x_max)

    all_axes[1].plot(time_steps, pos_z, '-bo')
    all_axes[1].grid(True)
    all_axes[1].set_xlabel("time step")
    all_axes[1].set_ylabel("position z [km]")
    all_axes[1].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[1].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[1].set_xlim(1, start_step_max)
    all_axes[1].set_ylim(pos_z_min, pos_z_max)

    all_axes[2].plot(time_steps, velo_x, '-bo')
    all_axes[2].grid(True)
    all_axes[2].set_xlabel("time step")
    all_axes[2].set_ylabel("velocity x [mm/yr]")
    all_axes[2].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[2].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[2].set_xlim(1, start_step_max)
    all_axes[2].set_ylim(velo_x_min, velo_x_max)

    all_axes[3].plot(time_steps, velo_z, '-bo')
    all_axes[3].grid(True)
    all_axes[3].set_xlabel("time step")
    all_axes[3].set_ylabel("velocity z [mm/yr]")
    all_axes[3].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[3].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[3].set_xlim(1, start_step_max)
    all_axes[3].set_ylim(velo_z_min, velo_z_max)

    all_axes[4].plot(time_steps, temperature, '-bo')
    all_axes[4].grid(True)
    all_axes[4].set_xlabel("time step")
    all_axes[4].set_ylabel("temperature [°C]")
    all_axes[4].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[4].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[4].set_xlim(1, start_step_max)
    all_axes[4].set_ylim(temperature_min, temperature_max)

    all_axes[5].plot(time_steps, time_values, '-bo')
    all_axes[5].grid(True)
    all_axes[5].set_xlabel("time step")
    all_axes[5].set_ylabel("time values [MYrs]")
    all_axes[5].yaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[5].xaxis.set_major_formatter(mti.FormatStrFormatter("%.2f"))
    all_axes[5].set_xlim(1, start_step_max)
    all_axes[5].set_ylim(time_values_min, time_values_max)

    figure_file = "{}_{}".format(surface_node, velocity_file_name.replace(".dat", ".png"))

    fig.savefig(figure_file, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("usage:")
        print("{} node_id".format(sys.argv[0]))
        sys.exit(1)
    else:
        surface_node = int(sys.argv[1])

        file_index = 1

        while True:
            velocity_file_name = "velocity_info_{0:04d}.bin".format(file_index)
            temperature_file_name = "time_temperature_history_{0:04d}.bin".format(file_index)

            if os.path.exists(velocity_file_name) and os.path.exists(temperature_file_name):
                process_file(velocity_file_name, temperature_file_name, surface_node)
                file_index += 1
            else:
                print("Finished!")
                break
