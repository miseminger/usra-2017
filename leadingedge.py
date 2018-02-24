#find neighbours and clustering for all time frames; make a .csv
#to this .csv, add the relevant cellprofiler data and velocity data
#for now, assume the cells are all moving in the x direction
#find the cell furthest to the right
#take all cells with x values from the furthest right to distance d to the left of it
#have the option to take several (n) strips of cells: total cell group width/n
#for each strip, labeled 1,2,3,4,...,n, find the neighbours 
#I'm not sure if this will work with Matlab or not!

#!/usr/bin/python2

"""
Author: Cindy Tan
Last Modified Date: 18 July 2017

Processes a MATLAB simulation csv file and determines cell neighbours at all timeframes
Arguments:
    -C CellProfiler: add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r: path to green and red CSV files (at least one file must be provided): use files such as /home/labmember/RoskelleyLab/EXP06-MeONGFP+EpENmCherry-40-60/15_CellProfiler/Output-Green/AllMyExpt_Cells_Green.csv
    -m or -p: radius in microns or pixels
    -c cluster cells: add flag if cells should be clustered and results shown on graph
    -o overlay: add flag if graph should be overlaid on a background image (provide image path as argument)
    -d print distance: add flag if graph should be annotated with distances
Outputs:
    - csv file of total, green, and red neighbours of all cells over all timeframes
"""

from datetime import datetime
start_time = datetime.now()
time_mark = start_time

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
        ' must be specified\n  > all other options are not required'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='green_file', help='path to green CSV file')
parser.add_option('-r', type='string', dest='red_file', help='path to red CSV file')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-C', action='store_true', dest='cp', help='use CellProfiler file instead of MATLAB')

(options, args) = parser.parse_args()
error = 0
if options.cp:
    microns_per_pixel = 0.8
    start_count = 0
    data_cols = ['Metadata_FrameNumber', 'Location_Center_X', 'Location_Center_Y','TrackObjects_Displacement_15','TrackObjects_DistanceTraveled_15','TrackObjects_TrajectoryX_15','TrackObjects_TrajectoryY_15','AreaShape_Area']
if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m)
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if not(options.green_file and options.red_file):
    error = 1
    print 'Error: no file input, use -g and -r'
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

## EXTRACT DATA
def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True, usecols=data_cols)
    return da.view((float, len(da.dtype.names)))

try:
    green, red = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')
max_time = int(green[-1, 0])


# Data extraction
def parse_csv(file):
    df = pd.read_csv(file)
    if options.cp:
        frame_col = 'Metadata_FrameNumber'
        data_cols = ['ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
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
    return df

if options.green_file and options.red_file:
    green_data, red_data = parse_csv(options.green_file), parse_csv(options.red_file)
    split = green_data.shape[0]  #split is the number of rows in green_data
    data = pd.concat([green_data, red_data]) #stack green_data on top of red_data and call this 'data'

data = data.set_index(np.arange(data.shape[0])) #set the index from 0 to one less than the number of rows in green_data

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
    return cell_type(index, 'str') + str(int(data.iloc[index]['Cell_ID']))

# Collect data for histogram and clustering
degrees = pd.DataFrame(0, index=np.arange(data.shape[0]), columns=['Cell_ID', 'Total', 'Green', 'Red']) #make an empty dataframe for storing everything
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

        if distance < float(radius):
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

        
# Script complete text
print '\n.csv file created: ' + csv_name
runtime = str(datetime.now() - start_time).split(':')
if float(runtime[0]) == 0:
    if float(runtime[1]) == 0:
        runtime = str(round(float(runtime[2]), 3)) + 's'
    else:
        runtime = runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
else:
    runtime = runtime[0] + 'h ' + runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
print 'Runtime: ' + runtime + '\n'        
