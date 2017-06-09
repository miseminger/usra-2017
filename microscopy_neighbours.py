#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 9 June 2017

Processes a MATLAB simulation csv file and determines cell neighbours at a specific timeframe
Arguments:
    -c CellProfiler - add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r - path to green and red CSV files
    -t timeframe - time at which cells are analyzed, specified as frame number
    -m or -p - radius in microns or pixels
    -d print distance - add flag if graph should be annotated with distances
Output: determines cell neighbours from radius, produces a graph of cells and neighbours and histogram of # of cells vs.
# of neighbours, saves both images
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from optparse import OptionParser
from graphviz import Graph

# Options parsing
parser = OptionParser()
parser.add_option('-c', '--cellprofiler', action='store_true', dest='cp', help='use CellProfiler file instead of MATLAB')
parser.add_option('-g', '--greenfile', type='string', dest='green_file', help='path to green CSV file')
parser.add_option('-r', '--redfile', type='string', dest='red_file', help='path to red CSV file')
parser.add_option('-t', '--time', type='int', dest='time', help='timeframe at which cells are analyzed')
parser.add_option('-m', '--radius_m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', '--radius_p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-d', '--printdist', action='store_true', dest='printdist', help='print distance on edges')

(options, args) = parser.parse_args()
if options.cp:
    microns_per_pixel = 0.8
else:
    microns_per_pixel = 0.4
if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m)
printdist = 0
if options.printdist:
    printdist = options.printdist

# Data extraction
def parse_csv(file):
    df = pd.read_csv(file)
    if options.cp:
        frame_col = 'Metadata_FrameNumber'
        data_cols = ['ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
    else:
        frame_col = 'Frame'
        data_cols = ['CellID', 'CentroidX', 'CentroidY']
    indices = np.where(df[frame_col] == options.time)[0]
    df = df[indices[0]:indices[-1] + 1][data_cols]
    df.columns = ['Cell', 'X', 'Y']
    return df

green_data, red_data = parse_csv(options.green_file), parse_csv(options.red_file)
split = green_data.shape[0]
data = pd.concat([green_data, red_data])

count = np.zeros(data.shape[0])

# Graphviz setup
scale = 0.025
def name():
    return 'neighbours_r' + str(int(radius)) + '_t' + str(options.time)

redf = '#FF9090'
greenf = '#75CF90'
redl = '#FF3535'
greenl = '#019A2F'

g = Graph(name=name(), format='png', engine='neato')

g.attr('node', pin='true', shape='circle', style='filled', width='0.5', fixedsize='true', fontname='courier')
g.attr('edge', fontsize='10.0', penwidth='2', fontname='courier')

# Plot graph and collect histogram data
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
    position = str(x * scale) + ',' + str(y * scale)

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

count = np.split(count, [split])

# Plot histogram
plt.hist([count[0], count[1]], np.arange(min(np.concatenate((count[0], count[1]))) - 0.5,
                                         max(np.concatenate((count[0], count[1]))) + 1.5, 1),
         lw=0, color=[greenf, redf], label=['Green', 'Red'])
plt.legend()
plt.xticks(np.arange(0, max(np.concatenate((count[0], count[1]))) + 1, 1))
plt.xlim([min(np.concatenate((count[0], count[1]))) - 0.5, max(np.concatenate((count[0], count[1]))) + 0.5])
plt.xlabel('# of neighbours')
plt.ylabel('# of cells')
plt.axes().yaxis.grid()
plt.savefig(name() + '.png')

g.render()
