#!/usr/bin/python2

"""
Author: Cindy Tan
Last Modified Date: 4 Aug 2017

Processes a MATLAB simulation csv file and determines cell neighbours at all frames
Arguments:
    -C CellProfiler: add flag if data is from CellProfiler (as opposed to MATLAB, default)
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    CSV file with one row for each timeframe with frame number, number of green cells, number of red cells, green-green neighbours, red-red neighbours, green-red neighbours
    Plot of number of green cells, number of red cells, green-green neighbours, red-red neighbours, green-red neighbours versus frame number
"""
import os

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

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
if options.cp:
    microns_per_pixel = 0.8
    start_count = 0
    data_cols = ['Metadata_FrameNumber', 'Location_Center_X', 'Location_Center_Y']
else:
    microns_per_pixel = 0.8
    start_count = 1
    data_cols = ['Frame', 'CentroidX', 'CentroidY']
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
    da = np.genfromtxt(file, delimiter=',', names=True, usecols=data_cols)  #make an array out of the .csv, called da
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    green, red = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')
max_time = int(green[-1, 0])  #max_time is the last frame number (last row, first column of green)

    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'
print 'Processing frames...'
print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')


#FIND NEIGHBOURS
def find_neighbours(primary, secondary):
    
    if primary is not secondary:
        ai = 5  #put red-green neighbours in column 5
    elif primary is green:
        ai = 3  #put green-green neighbours in column 3
    else:
        ai = 4  #put red-red neighbours in column 4
        
    time = 0
    while time <= primary[primary.shape[0],0]:  # while time <= last timeframe
            timeslist = primary[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
            firsttime = timeslist.index(time)  #index for first instance of frame number
            lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
        while i == time: #go through all the objects in the green array for a certain timeframe
            x, y = primary.iloc[i]['Location_Center_X'], primary.iloc[i]['Location_Center_Y']
   
            #now go through and find all the green neighbours of cell i in that same timeframe (these are called ni)
            if primary == secondary:
                ni_array = secondary[firsttime:i] + secondary[(i + 1):(lasttime + 1)]   #skip row i
            else:
                ni_array = secondary[firsttime:(lasttime + 1)]   #iterate over all objects
            for ni in ni_array: 
                nx, ny = secondary.iloc[ni]['Location_Center_X'], secondary.iloc[ni]['Location_Center_Y']
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                
                if distance < float(radius):
                    np_neighbours[i, ai] += 1  #increase the neighbour count in row i, column ai by one
            i += 1   
        time += 1

        
time = start_count
red_s = red
total_frames = max_time - start_count + 1 #total_frames is the number of frames remaining to process
np_neighbours = np.empty((total_frames, 6))  #why total_frames, not max_time?

#now loop through and find the neighbours for everything
while time < max_time:
    if time > 0 and time % 10 == 0:
        time_mark = mark()  #this is just for the timestamp
    np_neighbours[time, 0] = time  #populate the leftmost column of np_neighbours with the frame numbers as you go

    first = True
    while green[0, 0] == time: #start from the top with the first frame
        x, y = green[0, 1], green[0, 2]  #x and y are for the first cell
        green = np.delete(green, 0, 0) #remove the top row of green, ie. the first cell, so we can look at everything else
        
        if first:
            np_neighbours[time, 1], np_neighbours[time, 2] = find_neighbours(green, green, time) + 1, find_neighbours(green, red_s, time)
            first = False
        else:
            find_neighbours(green, green, time)
            find_neighbours(green, red_s, time)
    
    while red[0, 0] == time:
        x, y = red[0, 1], red[0, 2]
        red = np.delete(red, 0, 0)
        
        find_neighbours(red, red, time)
    
    red_s = np.delete(red_s, np.arange(np_neighbours[time, 2]), 0)
    time += 1


csv_name = 'neighbours_1' + '.csv'
count = 1
while os.path.isfile(csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    csv_name = 'neighbours_' + str(count) + '.csv'
np.savetxt(csv_name, np_neighbours, delimiter=',')

#plot graph of average neighbours over time
input_file = csv_name
frame_num, green_cells, red_cells, green_green_neighbours, red_red_neighbours, green_red_neighbours = np.loadtxt(input_file, delimiter=',', unpack=True, skiprows=0)

plt.plot(frame_num, green_cells, linestyle='-', label='green cells', color='g')
plt.plot(frame_num, red_cells, linestyle='-', label='red cells', color='r')
plt.plot(frame_num, green_green_neighbours, label='green-green neighbours', color='c')
plt.plot(frame_num, red_red_neighbours, label='red-red neighbours', color='m')
plt.plot(frame_num, green_red_neighbours, label='green-red neighbours', color='k')
plt.xlabel('frame #')
plt.ylabel('number')
plt.legend(loc='upper right',prop={'size':7})
plt.grid()
plt.savefig('neighbours_' + str(count) + '.png')

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print 'CSV produced: ' + csv_name
print 'Plot produced: ' + 'neighbours_' + str(count) + '.png'
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'


