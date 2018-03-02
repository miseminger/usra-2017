#!/usr/bin/python2

"""
Author: Cindy Tan
Last Modified Date: 4 Aug 2017

Processes a MATLAB simulation csv file and determines cell neighbours at all frames
Arguments:
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    CSV file with one row for each cell with frame number, object number, number of green-green neighbours, red-red neighbours, green-red neighbours
    """

import os

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import read_csv

from optparse import OptionParser

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) and number of timeframes (-t) are required\n  > exactly one radius parameter (-m/-p) must be specified\n  > -C is an optional flag'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-C', action='store_true', dest='cp',
                  help='use CellProfiler file instead of MATLAB')

(options, args) = parser.parse_args()
error = 0
#if options.cp:
microns_per_pixel = 0.8
start_count = 0
data_cols = ['Metadata_FrameNumber', 'ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
#else:
#    microns_per_pixel = 0.8
#    start_count = 1
#    data_cols = ['Frame', 'CentroidX', 'CentroidY']
if not(options.greenfile and options.redfile):
    error = 1
    print 'Error: insufficient file input, both files must be provided using -g and -r'

if options.radius_p:
    radius = float(options.radius_p)
elif options.radius_m:
    radius = float(options.radius_m) / microns_per_pixel
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if error:
    sys.exit('Use -h for options usage help')


## SETUP
print 'Setting up...'

# Time functions
def timestring(time):
    time = str(time).split(':')
    seconds = '{:6.2f}'.format(round(float(time[2]), 2)) + 's'
    minutes = '{:>4}'.format(time[1].lstrip('0') + 'm')
    hours = '{:>3}'.format(time[0] + 'h')
    if float(time[0]) == 0:
        if float(time[1]) == 0:
            time = seconds
        else:
            time = minutes + seconds
    else:
        time = hours + minutes + seconds
    return '{:>14}'.format(time)

def mark():
    now = dt.now()
    total = timestring(now - start_time)
    last_block = timestring(now - time_mark)
    print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
    return now


## EXTRACT DATA
def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True, usecols=data_cols, max_rows=5)  #make an array out of the .csv, called da
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    green, red = np_read_csv(options.greenfile), np_read_csv(options.redfile)
    #green[:,2:6] = green[:,2:6] * microns_per_pixel #convert all distances in green from pixels to microns
    #red[:,2:6] = red[:,2:6] * microns_per_pixel #convert all distances in green from pixels to microns
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')
max_time = int(green[-1, 0])  #max_time is the last frame number (last row, first column of green)

# Helper functions
def cell_type(index, format='int'):
    if index < green.shape[0]:
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

redf = '#FF9090'
greenf = '#75CF90'
redl = '#FF3535'
greenl = '#019A2F'

    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'
print 'Processing frames...'
print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')

total_frames = max_time - start_count + 1 #total_frames is the total number of frames :)

#FIND NEIGHBOURS
def find_neighbours(primary, secondary):

    if primary is not secondary:
        ai = 4  #put red-green (heterotypic) neighbours in column 5
    elif primary is green:
        ai = 2  #put green-green neighbours in column 3
    else:
        ai = 3  #put red-red neighbours in column 4
		
    time = start_count  #start at t=0 or t=1 for CellProfiler or Matlab
        
    while time < primary[(primary.shape[0] - 1),0] + 1:  # while time <= last timeframe
	print 'time'
	print time
        timeslist = primary[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        firsttime = timeslist.index(time)  #index for first instance of frame number
        lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
	#if time > 0 and time % 10 == 0: #this is just for the timestamp
            #time_mark = mark()  
	i = firsttime #start with the index of the first instance of the frame number
	for i in np.arange(firsttime,(lasttime + 1)):
	    print 'i:'
	    print i
            x, y = primary[i,2] * microns_per_pixel, primary[i,3] * microns_per_pixel
   
            #now go through and find all the green neighbours of cell i in that same timeframe (these are called ni)
	    if primary is secondary:
		ni_array = np.append(np.arange(firsttime,i), np.arange(i,(lasttime + 1)))
            else:
		ni_array = np.arange(firsttime,(lasttime + 1))
            for ni in ni_array: 
                nx, ny = secondary[ni,2] * microns_per_pixel, secondary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                
                if distance < float(radius):
                    np_neighbours[i, ai] += 1  #increase the neighbour count in row i, column ai by one
        time += 1

#now loop through and find the neighbours for everything
#set the matrix to put red neighbour information in
np_neighbours = np.zeros(((red.shape[0]), 5))  #first two columns are for frame number and cell ID, last 3 are for neighbours
np_neighbours[:,0] = red[:,0]  #first row of np_neighbours is Metadata_FrameNumber (frame #)
np_neighbours[:,1] = red[:,1]  #second row of np_neighbours is ObjectNumber (cell ID)
find_neighbours(red, red)
print 'find_neighbours(red, red)...'
print np_neighbours

find_neighbours(red, green)
print 'find_neighbours(red, green)...'
print np_neighbours
np_neighbours = np.concatenate((np_neighbours, np.delete(red,[0,1],1)), axis=1)   #add the rest of the datacols to np_neighbours before saving it:
np_neighbours_red = np_neighbours

#make the green and red np_neighbours called something different
#set the matrix to put green neighbour information in
np_neighbours = np.zeros(((green.shape[0]), 5))  #first two columns are for frame number and cell ID, last 3 are for neighbours
np_neighbours[:,0] = green[:,0]  #first row of np_neighbours is Metadata_FrameNumber (frame #)
np_neighbours[:,1] = green[:,1]  #second row of np_neighbours is ObjectNumber (cell ID)
find_neighbours(green, green)
find_neighbours(green, red)    
np_neighbours = np.concatenate((np_neighbours, np.delete(green,[0,1],1)), axis=1)   #add the rest of the datacols to np_neighbours before saving it:
np_neighbours_green = np_neighbours

#refine cell IDs to include color tags
#def cell_name(index, color):
#    return cell_type(index, 'str') + str(int(color[index,1]))

#i = 0
#while i < red.shape[0]:
#	np_neighbours_red[i,1] = 'RED_' + str(i+1)
#	np_neighbours_green[i,1] = 'GREEN_' + str(i+1)
#	i = i + 1

#combine red and green neighbours: red on top, green below
np_neighbours_merged = np.concatenate((np_neighbours_red, np_neighbours_green), axis=0)

#add column labels
columnlabels = ['Metadata_FrameNumber', 'ObjectNumber', 'Green-Green Neighbours', 'Red-Red Neighbours', 'Red-Green Neighbours','Location_Center_X (px)', 'Location_Center_Y (px)']
#np_neighbours_merged = pd.DataFrame(np_neighbours_merged, columns=columnlabels)
np_neighbours_green = pd.DataFrame(np_neighbours_green, columns=columnlabels)
np_neighbours_red = pd.DataFrame(np_neighbours_red, columns=columnlabels)

green_csv_name = 'neighbours_1' + '_green.csv'
red_csv_name = 'neighbours_1' + '_red.csv'
count = 1
while os.path.isfile(red_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    red_csv_name = 'neighbours_' + str(count) + '.csv'
count = 1
while os.path.isfile(green_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    green_csv_name = 'neighbours_' + str(count) + '.csv'

np_neighbours_red.to_csv(red_csv_name, index=False, header=True, sep=' ')
np_neighbours_green.to_csv(green_csv_name, index=False, header=True, sep=' ')

#np.savetxt(csv_name, np_neighbours_merged, delimiter=',')

#add column labels to .csv: this way messes up the order
#df = read_csv(csv_name)
#df.columns = ['Metadata_FrameNumber', 'ObjectNumber', 'Green-Green Neighbours', 'Red-Red Neighbours', 'Red-Green Neighbours']
#df.to_csv(csv_name)

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print 'CSV produced: ' + csv_name
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'


