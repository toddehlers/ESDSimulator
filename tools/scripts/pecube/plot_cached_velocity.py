#!/usr/bin/env python3

import os.path
import sys
import glob
import re
import struct
import math

# import numpy

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

all_ages_files = glob.glob("cached_*.dat")

def plot_data(input_file, title):
    x_pos = []
    z_pos = []
    vx_values = []
    vz_values = []

    counter = 0

    (file_x, file_y, file_z) = struct.unpack("iii", input_file.read(4 * 3))

    for x in range(0, file_x):
        for y in range(0, file_y):
            for z in range(0, file_z):
                (px, py, pz, vx, vy, vz) = struct.unpack("dddddd", input_file.read(8 * 6))

                # if px >= 110.0 and px <= 130.0:
                #     if pz >= 100.0 and pz <= 120.0:
                #         x_pos.append(px)
                #         z_pos.append(pz)
                #         vx_values.append(vx)
                #         vz_values.append(vz)

                counter = counter + 1
                if counter > 100:
                    x_pos.append(px)
                    z_pos.append(pz)
                    vx_values.append(vx)
                    vz_values.append(vz)
                    counter = 0

            ignore = input_file.read(8) # ignore surface


    print("min vx: {}, max vx: {}".format(min(vx_values), max(vx_values)))
    print("min vz: {}, max vz: {}".format(min(vz_values), max(vz_values)))

    plt.figure(1).clear()
    plt.figure(1).subplots_adjust(top=0.90, bottom=0.10, left=0.10, right=0.90)

    ax = plt.subplot(111)

    plt.xlabel("x pos [km]")
    plt.ylabel("z pos [km]")

    plt.quiver(x_pos, z_pos, vx_values, vz_values, scale=200.0)

    ax.set_title("input file: {}".format(title))

    plt.savefig("{}.png".format(title), dpi=400)

    # with open("velo_field_{}.txt".format(title), "w") as f:
    #     for (px, pz, vx, vz) in zip(x_pos, z_pos, vx_values, vz_values):
    #         f.write("{} {} {} {}\n".format(px, pz, vx, vz))

def process_all_files(all_ages_files):
    for filename in all_ages_files:
        with open(filename, "rb") as input_file:
            print("open file: '{}'".format(filename))
            output_filename = os.path.splitext(filename)[0]
            plot_data(input_file, output_filename)
            # break


process_all_files(all_ages_files)
