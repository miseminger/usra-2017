#!/usr/bin/python2

"""
Author: Madeline Iseminger
Last Modified Date: 31 Mar 2018

Processes CellProfiler or MATLAB .csv files and determines homotypic and heterotypic cell neighbour IDs and positions for each cell at all frames
use .csv files of the type: AllMyExpt_MyCells_Red.csv, AllMyExpt_MyCells_green.csv
Arguments:
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    Two CSV files (one for red, one for green) with one row for each cell with frame number, object number, number of green-green neighbours, red-red neighbours, green-red neighbours, and the x and y positions of each of these
    columns of final CSV: 'Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)'
    """

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import os
import sys
import math
import numpy as np
import pandas as pd
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
    data_cols = ['Metadata_FrameNumber', 'ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
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
    da = np.genfromtxt(file, delimiter=',', names=True, usecols=data_cols)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    green, red = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')
max_time = int(green[-1, 0])
    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'

total_frames = max_time - start_count + 1 #total_frames is the total number of frames :)

#FIND NEIGHBOURS
def find_neighbours(primary, secondary):

    print 'Processing frames...'
    print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')
        
#    time_mark = dt.now()
    
    def timecollector(time):
        time_mark = dt.now()

        if time > 0 and time % 10 == 0:
            now = dt.now()
            total = timestring(now - start_time)
            last_block = timestring(now - time_mark)
            print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
            time_mark =  dt.now()

        timeslist = primary[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        firsttime = timeslist.index(time)  #index for first instance of frame number
        lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
      
        def makecellblock(i):
    
            homo_neighbour_ids = []
            homo_neighbour_x = []
            homo_neighbour_y = []
    
            hetero_neighbour_ids = []
            hetero_neighbour_x = []
            hetero_neighbour_y = []
    
            x, y = primary[i,2] * microns_per_pixel, primary[i,3] * microns_per_pixel
   
        #now go through and find all the green neighbours of cell i in that same timeframe (these are called ni)
        #find homotypic neighbours
            homo_ni_array = np.arange(firsttime,(lasttime + 1))

            for ni in homo_ni_array: 
                nx, ny = primary[ni,2] * microns_per_pixel, primary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                    
                if distance < float(radius) and ni != i:
                        ni_id = primary[ni,1]
                        homo_neighbour_ids.append(ni_id)
                        homo_neighbour_x.append(primary[ni,2])
                        homo_neighbour_y.append(primary[ni,3])
            
            if len(homo_neighbour_ids) == 0:
                homo_neighbour_ids.append(0)
                homo_neighbour_x.append(0)
                homo_neighbour_y.append(0)

        #find heterotypic neighbour ids
            secondary_timeslist = secondary[:,0].tolist()
            ni_firsttime = secondary_timeslist.index(time)  #index for first instance of frame number
            ni_lasttime = len(secondary_timeslist) - secondary_timeslist[::-1].index(time) - 1 #index for last instance of frame number
            hetero_ni_array = np.arange(ni_firsttime,(ni_lasttime + 1))

            for ni in hetero_ni_array: 
                nx, ny = secondary[ni,2] * microns_per_pixel, secondary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                    
                if distance < float(radius):
                        ni_id = secondary[ni,1]
                        hetero_neighbour_ids.append(ni_id)
                        hetero_neighbour_x.append(secondary[ni,2])
                        hetero_neighbour_y.append(secondary[ni,3])
                        
            if len(hetero_neighbour_ids) == 0:
                hetero_neighbour_ids.append(0)
                hetero_neighbour_x.append(0)
                hetero_neighbour_y.append(0)

            homo_neighbour_ids = np.transpose(np.asarray(homo_neighbour_ids))
            homo_neighbour_ids.shape = (homo_neighbour_ids.shape[0],1)
            homo_neighbour_x = np.transpose(np.asarray(homo_neighbour_x))
            homo_neighbour_x.shape = (homo_neighbour_x.shape[0],1)
            homo_neighbour_y = np.transpose(np.asarray(homo_neighbour_y))
            homo_neighbour_y.shape = (homo_neighbour_y.shape[0],1)
    
            hetero_neighbour_ids = np.transpose(np.asarray(hetero_neighbour_ids))
            hetero_neighbour_ids.shape = (hetero_neighbour_ids.shape[0],1)
            hetero_neighbour_x = np.transpose(np.asarray(hetero_neighbour_x))
            hetero_neighbour_x.shape = (hetero_neighbour_x.shape[0],1)
            hetero_neighbour_y = np.transpose(np.asarray(hetero_neighbour_y))
            hetero_neighbour_y.shape = (hetero_neighbour_y.shape[0],1)
    
            #now make a block for the first cell ID, primary[i,1]    
            cellblocklength = max(homo_neighbour_ids.shape[0],hetero_neighbour_ids.shape[0])
    
            #make the arrays the same length by adding zeros on the end of the shorter one
            difference = np.abs(homo_neighbour_ids.shape[0] - hetero_neighbour_ids.shape[0])
            addzeros = np.zeros((difference,1))
            if homo_neighbour_ids.shape[0] < hetero_neighbour_ids.shape[0]:
                homo_neighbour_ids = np.concatenate((homo_neighbour_ids, addzeros))
                homo_neighbour_x = np.concatenate((homo_neighbour_x, addzeros))
                homo_neighbour_y = np.concatenate((homo_neighbour_y, addzeros))
            else:
                hetero_neighbour_ids = np.concatenate((hetero_neighbour_ids, addzeros))
                hetero_neighbour_x = np.concatenate((hetero_neighbour_x, addzeros))
                hetero_neighbour_y = np.concatenate((hetero_neighbour_y, addzeros))
            
            cellblock = np.zeros((cellblocklength,2))
            cellblock[:,0] = time  #first column is time
            cellblock[:,1] = primary[i,1]  #second column is cell ID
    
            #after you make all the timeblocks the right length, append the neighbour ids!
            cellblock = np.concatenate((cellblock, homo_neighbour_ids), axis=1)  #third column is homotypic neighbour IDs
            cellblock = np.concatenate((cellblock, hetero_neighbour_ids), axis=1)  #third column is homotypic neighbour IDs
            
            mainpos = np.zeros((cellblocklength,2))
            mainpos[:,0] = x  #first column is x position of main cell
            mainpos[:,1] = y  #second column is y position of main cell
    
            cellblock = np.concatenate((cellblock, mainpos), axis=1)  #fifth column is x position of main cell, sixth is y position
            cellblock = np.concatenate((cellblock, homo_neighbour_x), axis=1)  #seventh column is x position of homo neighbour
            cellblock = np.concatenate((cellblock, homo_neighbour_y), axis=1)  #eighth column is y position of homo neighbour
            cellblock = np.concatenate((cellblock, hetero_neighbour_x), axis=1)  #ninth column is x position of hetero neighbour
            cellblock = np.concatenate((cellblock, hetero_neighbour_y), axis=1)  #tenth column is y position of hetero neighbour
    
            return cellblock

        allcellblocks = np.concatenate([makecellblock(x) for x in np.arange(firsttime,(lasttime + 1))])
        return allcellblocks

    alltimes = np.concatenate([timecollector(x) for x in np.arange(start_count,(max_time + 1))])
    return alltimes


#now loop through and find the neighbours for everything

#define column labels
columnlabels = ['Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']

print "Finding neighbours of red cells..."
np_neighbours_red = find_neighbours(red, green)
np_neighbours_red = pd.DataFrame(np_neighbours_red, columns=columnlabels)

print "Finding neighbours of green cells..."
np_neighbours_green = find_neighbours(green, red)    
np_neighbours_green = pd.DataFrame(np_neighbours_green, columns=columnlabels)

green_csv_name = 'allneighbours_1' + '_green.csv'
red_csv_name = 'allneighbours_1' + '_red.csv'

count = 1
while os.path.isfile(red_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    red_csv_name = 'allneighbours_' + str(count) + '_red.csv'

count = 1
while os.path.isfile(green_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    green_csv_name = 'allneighbours_' + str(count) + '_green.csv'

np_neighbours_red.to_csv(red_csv_name, index=False, header=True, sep=',')
np_neighbours_green.to_csv(green_csv_name, index=False, header=True, sep=',')

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print 'CSVs produced: ' + green_csv_name + ' and ' + red_csv_name
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'

