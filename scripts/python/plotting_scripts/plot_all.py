import argparse
import matplotlib.pyplot as plt
import os
from collections import defaultdict

from scripts.python.plotting_scripts import plot_file_sizes as pfs
from scripts.python.plotting_scripts import plot_columns as pc

BLOCK_SIZES = [10000, 15000, 20000, 'map']

def parse_args():
    parser = argparse.ArgumentParser(description='Plot file sizes')
    parser.add_argument('--decomp', type=str, required=True,
                        help='dir with decompression times')
    parser.add_argument('--colors', type=str, required=True,
                        help='path to file which lists for each codec')
    parser.add_argument('--out', type=str, required=True,
                        help='Output directory for plots')
    return parser.parse_args()

def read_colors(colors_file):
    # codec, color
    colors = {}
    with open(colors_file, 'r') as f:
        for line in f:
            L = line.strip().split(',')
            colors[L[0]] = L[1]
    return colors


