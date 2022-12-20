#!/usr/bin/env python3

import os.path
import sys
import glob
import re
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

all_temperature_files = glob.glob("Temps_tec*.dat")
min_max_file = "min_max_values_temp.txt"

match_nodes = re.compile("(\s*)n(\s*)=(\s*)(\d+)(\s*),")
#  n=      519750 , e=      411584 et=brick,f=fepoint

mpl.rcParams["font.size"] = 6
mpl.rcParams["lines.markersize"] = 3.0
mpl.rcParams['axes.color_cycle'] = ["#0000ff", "#00aa00", "#aa0000", "#00aaaa", "#aa00aa", "#aaaa00", "#000000", "#8888ff"]



def process_file(input_file):
    num_of_entries = 0
    entry_counter = 0
    skip_counter = 0

    x_pos = []
    z_pos = []
    temperature = []

    for line in input_file:
        m = match_nodes.match(line)
        if m:
            num_of_entries = int(m.group(4))
            continue

        if (num_of_entries > 0):
            entries = [float(val) for val in line.split()]
            if entries[1] > 0.0:
                break
            if len(entries) >= 8:
                if skip_counter == 0:
                    x_pos.append(entries[0])
                    z_pos.append(entries[2])
                    temperature.append(entries[4])

                skip_counter += 1
                if skip_counter > 0:
                    skip_counter = 0

    # print("num_of_entries: {}".format(num_of_entries))

    x_min = min(x_pos)
    x_max = max(x_pos)
    y_min = min(z_pos)
    y_max = max(z_pos)

    return (x_pos, z_pos, temperature, x_min, x_max, y_min, y_max)

def plot_data(title, x_pos, z_pos, temperature, x_min, x_max, y_min, y_max):
    plt.figure(1).clear()
    plt.figure(1).subplots_adjust(top=0.96, bottom=0.06, left=0.06, right=0.98)

    ax = plt.subplot(111)

    plt.xlabel("x pos [km]")
    plt.ylabel("z pos [km]")
    plt.axis([x_min, x_max, y_min, y_max])

    sc = plt.scatter(x_pos, z_pos, c=temperature, edgecolors='none')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.xaxis.set_minor_locator(mti.MultipleLocator(5))
    ax.yaxis.set_minor_locator(mti.MultipleLocator(1))
    ax.set_title("input file: {}".format(title))

    plt.colorbar(sc)

    plt.savefig("{}.png".format(title), dpi=200)

def determin_min_max(all_temperature_files):
    global_x_min = sys.float_info.max
    global_x_max = sys.float_info.min
    global_y_min = sys.float_info.max
    global_y_max = sys.float_info.min

    for filename in all_temperature_files:
        with open(filename, "r") as input_file:
            print("open file: '{}'".format(filename))
            (x_pos, z_pos, temperature, x_min, x_max, y_min, y_max) = process_file(input_file)
            global_x_min = min(global_x_min, x_min)
            global_x_max = max(global_x_max, x_max)
            global_y_min = min(global_y_min, y_min)
            global_y_max = max(global_y_max, y_max)

    return (global_x_min, global_x_max, global_y_min, global_y_max)

def process_all_files(all_temperature_files, global_x_min, global_x_max,
    global_y_min, global_y_max):

    for filename in all_temperature_files:
        with open(filename, "r") as input_file:
            print("open file: '{}'".format(filename))
            (x_pos, z_pos, temperature, x_min, x_max, y_min, y_max) = process_file(input_file)
            plot_data(os.path.splitext(filename)[0],
                x_pos, z_pos, temperature, global_x_min,
                global_x_max, global_y_min, global_y_max)


if __name__ == "__main__":
    global_x_min = sys.float_info.max
    global_x_max = sys.float_info.min
    global_y_min = sys.float_info.max
    global_y_max = sys.float_info.min

    if os.path.exists(min_max_file):
        with open(min_max_file, "r") as f:
            (global_x_min, global_x_max, global_y_min, global_y_max) = [float(val) for val in f.readline().split()]
    else:
        (global_x_min, global_x_max, global_y_min, global_y_max) = determin_min_max(all_temperature_files)
        with open(min_max_file, "w") as f:
            f.write("{} {} {} {}\n".format(global_x_min, global_x_max, global_y_min, global_y_max))

    print("x_min: {}, x_max: {}, y_min: {}, y_max: {}".format(global_x_min, global_x_max, global_y_min, global_y_max))

    process_all_files(all_temperature_files, global_x_min, global_x_max, global_y_min, global_y_max)
