#!/usr/bin/env python3

import sys

def process_file(f):
    min_val = 100000.0
    max_val = -1000000.0

    for line in f:
        values = [float(v) for v in line.split()]
        min_val = min(min_val, min(values))
        max_val = max(max_val, max(values))

    print("min: {}, max: {}".format(min_val, max_val))


if __name__ == "__main__":
    args = sys.argv

    if len(args) != 2:
        print("usage: {} input_file".format(args[0]))
        sys.exit(0)

    filename = args[1]
    with open(filename, "r") as f:
        process_file(f)
