#!/usr/bin/python2

"""
Author: Cindy Tan
Last Modified Date: 20 July 2017

Processes a MATLAB simulation csv file and determines cell neighbours at all frames
Arguments:
    -C CellProfiler: add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Output:
    CSV file with one row for each timeframe with frame number, number of green cells, number of red cells, green-green neighbours, red-red neighbours, green-red neighbours
"""

from datetime import datetime
start_time = datetime.now()

import os
import sys
import math
import numpy as np
import pandas as pd

from optparse import OptionParser

# Options parsing
usage = '%prog options\n  > at least one file (-r/-g) must be specified\n  > exactly one radius parameter (-m/-p)' \
        ' must be specified\n  > time (-t) must be specified\n  > all other options are not required'
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
if error:
    sys.exit('Use -h for options usage help')

# Create output DataFrame

print 'Setting up...'
neighbours = pd.DataFrame(columns=['FrameNumber', 'GreenCount', 'RedCount', 'Green-Green', 'Red-Red', 'Green-Red'])
        
# Data extraction

if options.cp:
    frame_col = 'Metadata_FrameNumber'
    data_cols = ['ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
else:
    frame_col = 'Frame'
    data_cols = ['CellID', 'CentroidX', 'CentroidY']

def parse_csv(file):
    all_data = pd.read_csv(file)
    try:
        indices = np.where(all_data[frame_col] == 0)[0]
    except KeyError:
        sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
                 'Use -h for options usage help')
    return all_data

def extract_time(all_data, time):
    indices = np.where(all_data[frame_col] == time)[0]
    data = all_data[indices[0]:indices[-1] + 1][data_cols]
    data.columns = ['Cell_ID', 'X', 'Y']
    return data

# Body function to go in the loop, finds neighbours

def find_neighbours(data, split, time):
    
    def cell_type(index):
        if index < split:
            return 1
        else:
            return 2
            
    red_red = 0
    green_green = 0
    green_red = 0

    i = 0
    while i < data.shape[0]:
        x, y = data.iloc[i]['X'] * microns_per_pixel, data.iloc[i]['Y'] * microns_per_pixel

        ni = i + 1
        while ni < data.shape[0]:
            nx, ny = data.iloc[ni]['X'] * microns_per_pixel, data.iloc[ni]['Y'] * microns_per_pixel
            distance = math.sqrt((x - nx)**2 + (y - ny)**2)

            if distance < float(radius):
                
                if cell_type(i) == cell_type(ni):
                    if cell_type(i) == 1:
                        green_green += 1
                    else:
                        red_red += 1
                else:
                    green_red += 1
                    
            ni += 1
        i += 1

    neighbours.loc[len(neighbours)] = [time, split, data.shape[0] - split, green_green, red_red, green_red]

def status(time):
    runtime = str(datetime.now() - start_time).split(':')
    if float(runtime[0]) == 0:
        if float(runtime[1]) == 0:
            runtime = str(round(float(runtime[2]), 3)) + 's'
        else:
            runtime = runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
    else:
        runtime = runtime[0] + 'h ' + runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
    return '{:20}'.format(str(time)) + str(runtime)

# Loops
if options.green_file and options.red_file:
    green_all, red_all = parse_csv(options.green_file), parse_csv(options.red_file)
    print 'Processing frames...'
    print '{:20}'.format('Frames processed') + 'Runtime'
    time = 0
    while not np.where(green_all[frame_col] == time)[0].size == 0:
        if time % 10 == 0:
            print status(time)
        green, red = extract_time(green_all, time), extract_time(red_all, time)
        split = green.shape[0]
        df = pd.concat([green, red])
        find_neighbours(df, split, time)
        time += 1
else:
    if options.green_file:
        green_all = parse_csv(options.green_file)
        
        time = 0
        while not np.where(green_all[frame_col] == time)[0].size == 0:
            if time % 10 == 0:
                print status(time)
            df = extract_time(green_all, time)
            split = data.shape[0]
            find_neighbours(df, split, time)
            time += 1
        
    else:
        red_all = parse_csv(options.red_file)
        split = 0
        
        time = 0
        while not np.where(red_all[frame_col] == time)[0].size == 0:
            if time % 10 == 0:
                print status(time)
            df = extract_time(red_all, time)
            find_neighbours(df, split, time)
            time += 1

print '\n', neighbours, '\n'
neighbours.to_csv('neighbours.csv', index=False)

runtime = str(datetime.now() - start_time).split(':')
if float(runtime[0]) == 0:
    if float(runtime[1]) == 0:
        runtime = str(round(float(runtime[2]), 3)) + 's'
    else:
        runtime = runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
else:
    runtime = runtime[0] + 'h ' + runtime[1] + 'm ' + str(round(float(runtime[2]), 3)) + 's'
print 'Total runtime: ' + runtime + '\n'