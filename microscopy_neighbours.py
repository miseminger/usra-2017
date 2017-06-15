#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 15 June 2017

Processes a MATLAB simulation csv file and determines cell neighbours at a specific timeframe
Arguments:
    -c CellProfiler: add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -t timeframe: time at which cells are analyzed, specified as frame number
    -m or -p: radius in microns or pixels
    -d print distance: add flag if graph should be annotated with distances
    -o overlay: add flag if graph should be overlaid on a background image (provide image path as argument)
Output: determines cell neighbours from radius, produces a graph of cells and neighbours and histogram of # of cells vs.
# of neighbours, saves both images
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PIL import Image
from optparse import OptionParser
from graphviz import Graph

# Options parsing
usage = '%prog options\n  > at least one file (-r/-g) must be specified\n  > exactly one radius parameter (-m/-p)' \
        ' must be specified\n  > time (-t) must be specified\n  > all other options are not required'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='green_file', help='path to green CSV file')
parser.add_option('-r', type='string', dest='red_file', help='path to red CSV file')
parser.add_option('-t', type='int', dest='time', help='timeframe at which cells are analyzed')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-c', action='store_true', dest='cp', help='use CellProfiler file instead of MATLAB')
parser.add_option('-d', action='store_true', dest='printdist', help='print distance on edges')
parser.add_option('-o', type='string', dest='overlay', help='use graph as overlay, argument is path to background image file')

(options, args) = parser.parse_args()
error = 0
if options.printdist:
    printdist = options.printdist
else:
    printdist = 0
if options.cp:
    microns_per_pixel = 0.8
else:
    microns_per_pixel = 0.8
if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m)
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if not(options.green_file or options.red_file):
    error = 1
    print 'Error: no file input, use -g and/or -r'
if not(options.time):
    error = 1
    print 'Error: no timeframe specified, use -t'
if error:
    sys.exit('Use -h for options usage help')

# Store output information in a string
if options.green_file and options.red_file:
    info = '_r' + str(int(radius)) + '_t' + str(options.time)
else:
    if options.green_file:
        info = '_green_r' + str(int(radius)) + '_t' + str(options.time)
    else:
        info = '_red_r' + str(int(radius)) + '_t' + str(options.time)

# Build directory to contain outputs
directory = 'microscopy_neighbours' + info
directory_i = 0
while os.path.exists(directory):
    directory_i += 1
    directory = 'microscopy_neighbours' + info + '_' + str(directory_i)
os.makedirs(directory)

# Data extraction
def parse_csv(file):
    df = pd.read_csv(file)
    if options.cp:
        frame_col = 'Metadata_FrameNumber'
        data_cols = ['ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
    else:
        frame_col = 'Frame'
        data_cols = ['CellID', 'CentroidX', 'CentroidY']
    try:
        indices = np.where(df[frame_col] == options.time)[0]
    except KeyError:
        sys.exit('Error: could not parse CSV, use -c if file is for CellProfiler (MATLAB files are default)\n'
                 'Use -h for options usage help')
    df = df[indices[0]:indices[-1] + 1][data_cols]
    df.columns = ['Cell', 'X', 'Y']
    return df

if options.green_file and options.red_file:
    green_data, red_data = parse_csv(options.green_file), parse_csv(options.red_file)
    split = green_data.shape[0]
    data = pd.concat([green_data, red_data])
else:
    if options.green_file:
        data = parse_csv(options.green_file)
        split = data.shape[0]
    else:
        data = parse_csv(options.red_file)
        split = 0

count = np.zeros(data.shape[0])

# Graphviz setup
scale = 0.025

redf = '#FF9090'
greenf = '#75CF90'
redl = '#FF3535'
greenl = '#019A2F'

g = Graph(name='cellgraph' + info, format='png', engine='neato')

g.attr('node', pin='true', shape='circle', style='filled', width='0.5', fixedsize='true', fontname='courier')
g.attr('edge', fontsize='10.0', penwidth='2', fontname='courier')

# Create graph and collect histogram data
i = 0
while i < data.shape[0]:
    def cell_type(index, out='int'):
        if index < split:
            if out == 'str':
                return 'G'
            if out == 'int':
                return 1
            if out == 'hex':
                return greenf
        else:
            if out == 'str':
                return 'R'
            if out == 'int':
                return 2
            if out == 'hex':
                return redf

    def cell_name(index):
        return cell_type(index, 'str') + str(int(data.iloc[index]['Cell']))

    x, y = data.iloc[i]['X'] * microns_per_pixel, data.iloc[i]['Y'] * microns_per_pixel
    position = str(x * scale) + ',' + str(y * scale * -1)

    g.node(cell_name(i), pos=position, fillcolor=cell_type(i, 'hex'))

    ni = i + 1
    while ni < data.shape[0]:
        nx, ny = data.iloc[ni]['X'] * microns_per_pixel, data.iloc[ni]['Y'] * microns_per_pixel
        distance = math.sqrt((x - nx)**2 + (y - ny)**2)

        if distance < float(radius):
            # print cell_name(ni),
            if cell_type(i) == cell_type(ni):
                lcolor = cell_type(i, 'hex')
            else:
                lcolor = 'black'

            label = ''
            if printdist:
                label = str(round(distance, 2))
            g.edge(cell_name(i), cell_name(ni), color=lcolor, label=label)

            count[i] += 1
            count[ni] += 1

        ni += 1

    i += 1

# Generate graph image
g.render(directory + '/cellgraph' + info + '.gv')

# Plot and generate histogram image
count = np.split(count, [split])

plt.hist([count[0], count[1]], np.arange(min(np.concatenate((count[0], count[1]))) - 0.5,
                                         max(np.concatenate((count[0], count[1]))) + 1.5, 1),
         lw=0, color=[greenf, redf], label=['Green', 'Red'])
plt.legend()
plt.xticks(np.arange(0, max(np.concatenate((count[0], count[1]))) + 1, 1))
plt.xlim([min(np.concatenate((count[0], count[1]))) - 0.5, max(np.concatenate((count[0], count[1]))) + 0.5])
plt.xlabel('# of neighbours')
plt.ylabel('# of cells')
plt.axes().yaxis.grid()
plt.savefig(directory + '/neighbourhist' + info + '.png')

# Overlay image
if options.overlay:
    image = Image.open(options.overlay)
    graph = Image.open(directory + '/cellgraph' + info + '.gv.png')
    ratio = float(image.size[0]) / image.size[1]

    image = image.resize((int(round(graph.size[1] * ratio)), graph.size[1]))

    mask = Image.new('L', graph.size, color=0)

    bands = list(graph.split())
    if len(bands) == 3:
        mask = Image.new('L', graph.size, color=128)
        graph.putalpha(mask)
        bands = list(graph.split())
    else:
        bands[3] = bands[3].point(lambda x: x * 0.5)
    graph = Image.merge(graph.mode, bands)
    image.paste(graph, (0, 0), graph)

    image.save(directory + '/cellgraph_overlay' + info + '.png')
