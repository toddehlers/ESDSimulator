#!/usr/bin/env python3

import os.path
import sys
import glob
import re

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mti

all_ages_files = glob.glob("Ages_tec*.dat")
min_max_file = "min_max_values_position.txt"

match_variables = re.compile("VARIABLES(\s*)=(\s*)")
match_variable_name = re.compile("\"(.+?)\"\s*?")
match_nodes = re.compile("(\s*)n(\s*)=(\s*)(\d+)(\s*),")

mpl.rcParams["font.size"] = 6
mpl.rcParams["lines.markersize"] = 3.0
mpl.rcParams['axes.color_cycle'] = ["#0000ff", "#00aa00", "#aa0000", "#00aaaa", "#aa00aa", "#aaaa00", "#000000", "#8888ff"]

def process_file(input_file, offset):
    name_of_ages = []
    num_of_entries = 0
    num_of_ages = 0
    entry_counter = 0

    x_pos = []
    z_pos = []

    for line in input_file:
        m = match_variables.match(line)
        if m:
            name_of_ages = match_variable_name.findall(line)[(3 + offset):]
            num_of_ages = len(name_of_ages)
            continue
        m = match_nodes.match(line)
        if m:
            num_of_entries = int(m.group(4))
            continue

        if (num_of_entries > 0) and (num_of_ages > 0):
            entries = [float(val) for val in line.split()]
            if entries[1 + offset] > 0.0:
                break

            if len(entries) == num_of_ages + 3 + offset:
                x_pos.append(entries[offset])

                z_pos.append(entries[(2 + offset)])

                entry_counter += 1
                if entry_counter == num_of_entries:
                    break

    print("entries: {}".format(num_of_entries))

    x_min = min(x_pos)
    x_max = max(x_pos)
    z_min = min(z_pos)
    z_max = max(z_pos)

    return (x_pos, z_pos, name_of_ages, x_min, x_max, z_min, z_max)

def plot_data(title, x_pos, z_pos, name_of_ages, x_min, x_max, y_min, y_max):
    plt.figure(1).clear()
    plt.figure(1).subplots_adjust(top=0.96, bottom=0.06, left=0.06, right=0.98)

    ax = plt.subplot(111)

    plt.xlabel("x pos [km]")
    plt.ylabel("z pos [km]")
    plt.axis([x_min, x_max, y_min, y_max])

    plt.plot(x_pos, z_pos, "+")

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.xaxis.set_minor_locator(mti.MultipleLocator(5))
    #ax.yaxis.set_minor_locator(mti.MultipleLocator(1))
    ax.set_title("input file: {}".format(title))

    plt.savefig("Position_{}.png".format(title), dpi=200)

def determin_min_max(all_ages_files, offset):
    global_x_min = sys.float_info.max
    global_x_max = sys.float_info.min
    global_y_min = sys.float_info.max
    global_y_max = sys.float_info.min

    for filename in all_ages_files:
        with open(filename, "r") as input_file:
            print("open file: '{}'".format(filename))
            (x_pos, z_pos, name_of_ages, x_min, x_max, y_min, y_max) = process_file(input_file, offset)
            global_x_min = min(global_x_min, x_min)
            global_x_max = max(global_x_max, x_max)
            global_y_min = min(global_y_min, y_min / 1.1)
            global_y_max = max(global_y_max, y_max * 1.1)

    return (global_x_min, global_x_max, global_y_min, global_y_max)

def process_all_files(all_ages_files, global_x_min, global_x_max,
    global_y_min, global_y_max, offset):

    for filename in all_ages_files:
        with open(filename, "r") as input_file:
            print("open file: '{}'".format(filename))
            (x_pos, z_pos, name_of_ages, x_min, x_max, y_min, y_max) = process_file(input_file, offset)
            plot_data(os.path.splitext(filename)[0],
                x_pos, z_pos, name_of_ages, global_x_min,
                global_x_max, global_y_min, global_y_max)

if __name__ == "__main__":
    global_x_min = sys.float_info.max
    global_x_max = sys.float_info.min
    global_y_min = sys.float_info.max
    global_y_max = sys.float_info.min
    offset = 3

    if os.path.exists(min_max_file):
        with open(min_max_file, "r") as f:
            (global_x_min, global_x_max, global_y_min, global_y_max) = [float(val) for val in f.readline().split()]
    else:
        (global_x_min, global_x_max, global_y_min, global_y_max) = determin_min_max(all_ages_files, offset)
        with open(min_max_file, "w") as f:
            f.write("{} {} {} {}\n".format(global_x_min, global_x_max, global_y_min, global_y_max))

    print("x_min: {}, x_max: {}, y_min: {}, y_max: {}".format(global_x_min, global_x_max, global_y_min, global_y_max))

    process_all_files(all_ages_files, global_x_min, global_x_max, global_y_min, global_y_max, offset)
