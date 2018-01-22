#!/usr/bin/python2

"""
Author: Cindy Tan
Last Modified Date: 18 July 2017

Processes a MATLAB simulation csv file and determines cell neighbours at a specific timeframe
Arguments:
    -C CellProfiler: add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -t timeframe: time at which cells are analyzed, specified as frame number
    -m or -p: radius in microns or pixels
    -c cluster cells: add flag if cells should be clustered and results shown on graph
    -o overlay: add flag if graph should be overlaid on a background image (provide image path as argument)
    -d print distance: add flag if graph should be annotated with distances
Outputs:
    - histogram of total neighbours of red and green cells
    - csv file of total, green, and red neighbours of all cells
    - graph of cells with neighbour connections
        - graphviz (.gv) file for graph
        - can be clustered
    optional:
    - merged microscope image
    - overlay of graph on microscope image
"""

from datetime import datetime
start_time = datetime.now()

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PIL import Image, ImageMath
from optparse import OptionParser
from graphviz import Graph
from sklearn.cluster import DBSCAN

# Options parsing
usage = '%prog options\n  > at least one file (-r/-g) must be specified\n  > exactly one radius parameter (-m/-p)' \
        ' must be specified\n  > time (-t) must be specified\n  > all other options are not required'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='green_file', help='path to green CSV file')
parser.add_option('-r', type='string', dest='red_file', help='path to red CSV file')
parser.add_option('-t', type='int', dest='time', help='timeframe at which cells are analyzed')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-C', action='store_true', dest='cp', help='use CellProfiler file instead of MATLAB')
parser.add_option('-d', action='store_true', dest='printdist', help='print distance on edges')
parser.add_option('-c', action='store_true', dest='cluster', help='cluster cells')
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
if options.overlay:
    overlay = options.overlay.split(',')
    if len(overlay) > 2:
        error = 1
        print 'Error: too many image files for overlay'
if error:
    sys.exit('Use -h for options usage help')

# Store input information in a string
if options.cp:
    info = '_CP'
else:
    info = '_ML'
if not(options.green_file and options.red_file):
    if options.green_file:
        info += '_green'
    else:
        info += '_red'
info += '_r' + str(int(radius)) + '_t' + str(options.time)

# Build directory to contain outputs
directory = 'neighbours' + info
directory_i = 0
while os.path.exists(directory):
    directory_i += 1
    directory = 'neighbours' + info + '_' + str(directory_i)
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
        if not len(indices):
            max = df.loc[df.shape[0] - 1][frame_col]
            sys.exit('Error: Timeframe out of range\nMaximum timeframe: ' + str(max) + '\nTimeframe provided: ' + str(options.time))
    except KeyError:
        sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
                 'Use -h for options usage help')
    df = df[indices[0]:indices[-1] + 1][data_cols]
    df.columns = ['Cell_ID', 'X', 'Y']
    return df  #this is a DataFrame

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

data = data.set_index(np.arange(data.shape[0]))

# Helper functions
def cell_type(index, format='int'):
    if index < split:
        if format == 'fullstr':
            return 'Green'
        if format == 'str':
            return 'G'
        if format == 'int':
            return 1
        if format == 'hex':
            return greenf
    else:
        if format == 'fullstr':
            return 'Red'
        if format == 'str':
            return 'R'
        if format == 'int':
            return 2
        if format == 'hex':
            return redf

def cell_name(index):
    return cell_type(index, 'str') + str(int(data.iloc[index]['Cell_ID'])) #iloc fetches the Cell ID from row number (index)

# Collect data for histogram and clustering

# Code to determine intersections of lines drawn between centroids (from https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python)
from __future__ import division 

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

#make an empty chart with columns named as "columns", below, with the first column unlabeled but numbered 0 to n-1 and the rest filled with 0's
degrees = pd.DataFrame(0, index=np.arange(data.shape[0]), columns=['Cell_ID', 'TrueTotal', 'TrueGreen', 'TrueRed','Total', 'Green', 'Red'])  #  TrueTotal holds the number of true neighbours, TrueGreen the number of true neighbours that are green, and TrueRed the number of true neighbours that are red.

#fill the Cell_ID column with the Cell_ID numbers
degrees['Cell_ID'] = data['Cell_ID']

cell_positions = np.empty((2, data.shape[0]))
cell_neighbours = np.empty(data.shape[0], dtype=object)

i = 0
while i < data.shape[0]:

    cell_positions[0, i], cell_positions[1, i] = data.iloc[i]['X'] * microns_per_pixel, \
                                                 data.iloc[i]['Y'] * microns_per_pixel
    x, y = cell_positions[0, i], cell_positions[1, i]
    cell_neighbours[i] = []

    ni = i + 1
    while ni < data.shape[0]:
        nx, ny = data.iloc[ni]['X'] * microns_per_pixel, data.iloc[ni]['Y'] * microns_per_pixel
        distance = math.sqrt((x - nx)**2 + (y - ny)**2)
        
        #probably need to add another nesting thing--want to check against all other lines (check over all lines--one against all others)
        L1 = line([0,1], [2,3])
        L2 = line([2,3], [0,4])
        R = intersection(L1, L2)
        if R:
            print "Intersection detected:", R
        else:
            return False

        if distance < float(radius) and R == False:
            cell_neighbours[i].append(ni)

            degrees[cell_type(ni, 'fullstr')][i] += 1
            degrees[cell_type(i, 'fullstr')][ni] += 1

        ni += 1
    degrees['Total'][i] = degrees['Green'][i] + degrees['Red'][i]

    i += 1

degrees.to_csv(directory + '/cell_degrees' + info + '.csv', index=False)

# Clustering
if options.cluster:
    min_pts = int(math.ceil(np.mean(np.array(degrees['Total'])))) + 1

    db = DBSCAN(eps=radius, min_samples=min_pts).fit(np.transpose(cell_positions))

    num_colors = 12

    labels = db.labels_

    label_i = 0
    while label_i < len(db.labels_):
        label = labels[label_i]
        if label == -1:
            label = 0
        else:
            label = label % num_colors + 1
        labels[label_i] = label
        label_i += 1

# Generate graph image
scale = 0.025

redf = '#FF9090'
greenf = '#75CF90'
redl = '#FF3535'
greenl = '#019A2F'

g = Graph(name='cellgraph' + info, format='png', engine='neato')

g.attr('node', pin='true', shape='circle', width='0.5', fixedsize='true', style='filled', penwidth='3',
       colorscheme='set312', fontname='courier', fontsize='10.0')
g.attr('edge', fontsize='10.0', penwidth='2', fontname='courier')

i = 0
while i < data.shape[0]:
    x, y = cell_positions[0, i], cell_positions[1, i]
    position = str(x * scale) + ',' + str(y * scale * -1)

    if options.cluster:
        note = '\nC' + str(db.labels_[i])
        if labels[i] == 0:
            fill = 'white'
        else:
            fill = str(labels[i])
    else:
        note = ''
        fill = cell_type(i, 'hex')

    g.node(cell_name(i), label=cell_name(i) + note, pos=position, color=cell_type(i, 'hex'),
           fillcolor=fill)

    for ni in cell_neighbours[i]:
        if cell_type(i) == cell_type(ni):
            lcolor = cell_type(i, 'hex')
        else:
            lcolor = 'black'

        label = ''
        if printdist:
            label = str(round(distance, 2))

        g.edge(cell_name(i), cell_name(ni), color=lcolor, label=label)
    i += 1

g.render(directory + '/cellgraph' + info + '.gv')

# Plot and generate histogram image
degrees_split = np.split(np.array(degrees['Total']), [split])

plt.figure()
plt.hist([degrees_split[0], degrees_split[1]],
         np.arange(min(np.concatenate((degrees_split[0], degrees_split[1]))) - 0.5,
                   max(np.concatenate((degrees_split[0], degrees_split[1]))) + 1.5, 1),
         lw=0, color=[greenf, redf], label=['Green', 'Red'])
plt.legend()
plt.xticks(np.arange(0, max(np.concatenate((degrees_split[0], degrees_split[1]))) + 1, 1))
plt.xlim([min(np.concatenate((degrees_split[0], degrees_split[1]))) - 0.5,
          max(np.concatenate((degrees_split[0], degrees_split[1]))) + 0.5])
plt.xlabel('# of neighbours')
plt.ylabel('# of cells')
plt.axes().yaxis.grid()
plt.savefig(directory + '/neighbourhist' + info + '.png')

# Overlay image
if options.overlay:
    if len(overlay) == 1:
        image = Image.open(overlay[0])
    else:
        image0 = Image.open(overlay[0])
        image1 = Image.open(overlay[1])

        r0, g0, b0 = image0.split()[0], image0.split()[1], image0.split()[2]
        r1, g1, b1 = image1.split()[0], image1.split()[1], image1.split()[2]

        r, g, b = ImageMath.eval('convert(max(a, b), "L")', a=r0, b=r1), \
                  ImageMath.eval('convert(max(a, b), "L")', a=g0, b=g1), \
                  ImageMath.eval('convert(max(a, b), "L")', a=b0, b=b1)

        image = Image.merge('RGB', (r, g, b))
        image.save(directory + '/merged_image' + info + '.png')
    graph = Image.open(directory + '/cellgraph' + info + '.gv.png')
    ratio = float(image.size[0]) / image.size[1]

    image = image.resize((int(round(graph.size[1] * ratio)), graph.size[1]))

    bands = list(graph.split())
    if len(bands) == 3:
        mask = Image.new('L', graph.size, color=128)
        graph.putalpha(mask)
        bands = list(graph.split())
    else:
        bands[3] = bands[3].point(lambda x: x * 0.65)
    graph = Image.merge(graph.mode, bands)
    image.paste(graph, (0, 0), graph)

    image.save(directory + '/cellgraph_overlay' + info + '.png')

# Script complete text
print '\nOutput folder created: ' + directory
runtime = str(datetime.now() - start_time).split(':')
if float(runtime[0]) == 0:
    if float(runtime[1]) == 0:
        runtime = str(round(float(runtime[2]), 3)) + 's'
    else:
        runtime = runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
else:
    runtime = runtime[0] + 'h ' + runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
print 'Runtime: ' + runtime + '\n'
