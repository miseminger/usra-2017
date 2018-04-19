#!/usr/bin/python2

"""
Author: Madeline Iseminger
Last Modified Date: 31 Mar 2018

Implements a binary mask technique to find the leading edge of the cell mass in each frame, then increments by 35um along the leading edge 
and matches the position of each increment to the nearest cell found in the CellProfiler output .csv.  The IDs and positions of 
these leading edge cells are recorded in arrays (one for green cells, one for red) and the IDs and positions of the 
homotypic and heterotypic neighbours of the leading cells are then found and recorded in output .csv files.

Use .csv files of the type: AllMyExpt_MyCells_Red.csv, AllMyExpt_MyCells_green.csv, and phase images.

Arguments:
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    Two CSV files (one for red, one for green) with one row for each cell with frame number, object number, number of green-green neighbours, red-red neighbours, green-red neighbours, and the x and y positions of each of these
    columns of final CSV: 'Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)'
    """
    
#note:  for the CellProfiler data, (0,0) is in the lower left-hand corner of the image (http://forum.cellprofiler.org/t/xy-center-coordinates-for-objects-and-all-their-neighbors/627/3); Python indexes it as the upper lefthand corner.

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import os
import sys
import math
import numpy as np
import pandas as pd
from optparse import OptionParser

from scipy import ndimage
#from imageio import imread
#import imageio
from PIL import Image
from scipy.spatial import cKDTree
from skimage.filters.rank import entropy
from skimage.morphology import disk

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) and number of timeframes (-t) are required\n  > exactly one radius parameter (-m/-p) must be specified\n  > -C is an optional flag'
parser = OptionParser(usage=usage)
parser.add_option('-i', type='string', dest='imgfile', help='path to any of the ch00 phase images')
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
    
if options.imgfile:
    imgfile = options.imgfile    
if not(options.imgfile):
    error = 1
    print 'Error: insufficient file input, include the filepath to one of the phase images'
    

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

#max_time = int(green[-1, 0])    #max_time is the final frame number
max_time = 10
    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'

total_frames = max_time - start_count + 1 #total_frames is the total number of frames :)



def leadingarray_time():
    
    print 'Finding leading edge cells...'
    print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')
    
    def timeblock(time):
        
        time_mark = dt.now()
    
        if time > 0 and time % 10 == 0:
                now = dt.now()
                total = timestring(now - start_time)
                last_block = timestring(now - time_mark)
                print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
                time_mark =  dt.now()
                
        filename = imgfile[:(imgfile[0:-3].rindex('t') + 1)] + str("%03d" % (time,)) + imgfile[imgfile.rindex('_'):]
        
        greentimeslist = green[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        greenfirsttime = greentimeslist.index(time)  #index for first instance of frame number
        greenlasttime = len(greentimeslist) - greentimeslist[::-1].index(time) - 1 #index for last instance of frame number
        
        redtimeslist = red[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        redfirsttime = redtimeslist.index(time)  #index for first instance of frame number
        redlasttime = len(redtimeslist) - redtimeslist[::-1].index(time) - 1 #index for last instance of frame number
        
        leadingedge_green = []
        leadingedge_green_x = []
        leadingedge_green_y = []
        
        leadingedge_red = []
        leadingedge_red_x = []
        leadingedge_red_y = []
               
        #read the image into a numpy array called image_array
#        image_array = imageio.imread(filename)
#        image_array = np.asarray(image_array)

        im = Image.open(filename)
        image_array = np.array(im)
    
        max_y = image_array.shape[0]  #y is in the vertical direction
        
        #apply entropy mask    
        ent = entropy(image_array, disk(5))
        threshold = ent.mean()
        entmask = ent > threshold  #gives boolean array
        entmask = entmask * np.uint8(255) #cell area is mostly 255 (white), background is mostly 0 (black)
        
        #despeckle
        despeckled = ndimage.median_filter(entmask, size=50)  
        
        ##save image
        #imageio.imwrite('background.tif', image_array)
        #imageio.imwrite('overlay.tif', despeckled)
        
        #increment every 35-40um along the y axis
        stepsize = int(35. / microns_per_pixel)
        
#        if max_y % stepsize == 0:
#            extrabit = 1
#        else:
#            extrabit = max_y % stepsize
            
        for y in np.arange(0,max_y,stepsize):
            yslice = despeckled[y,:].tolist()  #list all the pixels for a given y, starting with y=0
            try:
                x = yslice[::1].index(0) - 1  #x is the white pixel right before the first black pixel when coming from the left
            except ValueError:
                continue   #if the wound has healed in that row, skip it and find the leading cell in the next row
                
            #may want a bit of code that counts up the rows in which the wound has healed for each frame
            
            # Shoe-horn existing data for entry into KDTree routines
            #  https://stackoverflow.com/questions/10818546/finding-index-of-nearest-point-in-numpy-arrays-of-x-and-y-coordinates?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
            points_list = [x,y]
        
            green_x_array = green[greenfirsttime:(greenlasttime + 1),2]
            green_y_array = green[greenfirsttime:(greenlasttime + 1),3]
        
            red_x_array = red[redfirsttime:(redlasttime + 1),2]
            red_y_array = red[redfirsttime:(redlasttime + 1),3]
            
            combined_green_x_y_arrays = np.dstack([green_y_array.ravel(),green_x_array.ravel()])[0]
            combined_red_x_y_arrays = np.dstack([red_y_array.ravel(),red_x_array.ravel()])[0]
            
            def do_kdtree(combined_x_y_arrays,points):
                mytree = cKDTree(combined_x_y_arrays)
                dist, indexes = mytree.query(points, k=1)  #find the closest point
                return dist, indexes
            
            greenresults = do_kdtree(combined_green_x_y_arrays,points_list)
            green_dist = greenresults[0]
            green_indexes = greenresults[1]
            
            redresults = do_kdtree(combined_red_x_y_arrays,points_list)
            red_dist = redresults[0]
            red_indexes = redresults[1]
            
            #if that position is within one radius of the white pixel, put that cell id in a list of leading edge cells (if it isn't already there)
            if red_dist < green_dist and red_dist < radius:
                leadingedge_red.append(red[red_indexes,1])
                leadingedge_red_x.append(red[red_indexes,2])
                leadingedge_red_y.append(red[red_indexes,3])
            elif green_dist < red_dist and green_dist < radius:
                leadingedge_green.append(green[green_indexes,1])
                leadingedge_green_x.append(green[green_indexes,2])
                leadingedge_green_y.append(green[green_indexes,3])
                
        greenarray = np.zeros((len(leadingedge_green),4))
        greenarray[:,0] = time  
        greenarray[:,1] = np.transpose(np.asarray(leadingedge_green))  #ids of green leading cells at time t
        greenarray[:,2] = np.transpose(np.asarray(leadingedge_green_x))   #x positions of green leading cells at time t
        greenarray[:,3] = 1200 - np.transpose(np.asarray(leadingedge_green_y)) #y positions of green leading cells at time t, transformed to match CP coordinate system
     
        redarray = np.zeros((len(leadingedge_red),4))
        redarray[:,0] = time  
        redarray[:,1] = np.transpose(np.asarray(leadingedge_red))  #ids of green leading cells at time t
        redarray[:,2] = np.transpose(np.asarray(leadingedge_red_x))   #x positions of red leading cells at time t
        redarray[:,3] = 1200 - np.transpose(np.asarray(leadingedge_red_y)) #y positions of red leading cells at time t, transformed to match CP coordinate system

        return greenarray, redarray

    allgreens = np.concatenate([timeblock(x)[0] for x in np.arange(0,(max_time + 1))])
    allreds = np.concatenate([timeblock(x)[1] for x in np.arange(0,(max_time + 1))])
    return allgreens, allreds

greens = leadingarray_time()[0]
reds = leadingarray_time()[1]


#FIND NEIGHBOURS
def find_neighbours(leadcolor, hetero_leadcolor, primary, secondary):

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

        #primary is the array of all red or green cells for all times
        timeslist = primary[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        firsttime = timeslist.index(time)  #index for first instance of frame number
        lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
        
        #leadcolor is the red or green array of leading edge cells
        leadcolortimeslist = leadcolor[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        leadcolorfirsttime = leadcolortimeslist.index(time)  #index for first instance of frame number
        leadcolorlasttime = len(leadcolortimeslist) - leadcolortimeslist[::-1].index(time) - 1 #index for last instance of frame number
        leaders_time = leadcolor[leadcolorfirsttime:(leadcolorlasttime + 1),1]
        
        hetero_leadcolortimeslist = hetero_leadcolor[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        hetero_leadcolorfirsttime = hetero_leadcolortimeslist.index(time)  #index for first instance of frame number
        hetero_leadcolorlasttime = len(hetero_leadcolortimeslist) - hetero_leadcolortimeslist[::-1].index(time) - 1 #index for last instance of frame number
        hetero_leaders_time = hetero_leadcolor[hetero_leadcolorfirsttime:(hetero_leadcolorlasttime + 1),1]
      
        def makecellblock(i):
            
            #first make a list of all the cell ids at the front
            #then go find all the neighbours for each of red and green
    
            homo_neighbour_ids = []
            homo_neighbour_x = []
            homo_neighbour_y = []
    
            hetero_neighbour_ids = []
            hetero_neighbour_x = []
            hetero_neighbour_y = []
    
            x, y = leadcolor[i,2] * microns_per_pixel, leadcolor[i,3] * microns_per_pixel
   
        #now go through and find all the green neighbours of cell i in that same timeframe (these are called ni)
        #find homotypic neighbours
            homo_ni_array = np.arange(firsttime,(lasttime + 1))

            for ni in homo_ni_array: 
                nx, ny = primary[ni,2] * microns_per_pixel, primary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                ni_id = primary[ni,1]
                
                #make sure to not include the leading cell itself as one of its neighbours!                        
                if (distance < float(radius)) and (ni_id != leadcolor[i,1]): #to exclude other leaders from the neighbours, change to "if distance < float(radius) and ni_id not in list(leaders_time):"
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
                ni_id = secondary[ni,1]
                nx, ny = secondary[ni,2] * microns_per_pixel, secondary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                
                if distance < float(radius): #to exclude other leaders from the neighbours, change to "if distance < float(radius) and ni_id not in list(list(hetero_leaders_time):"
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
    
            #now make a block for the first cell ID  
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
            cellblock[:,1] = leadcolor[i,1]  #second column is cell ID
    
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

        allcellblocks = np.concatenate([makecellblock(x) for x in np.arange(leadcolorfirsttime,(leadcolorlasttime + 1))])
        return allcellblocks

    alltimes = np.concatenate([timecollector(x) for x in np.arange(1,(max_time + 1))])
    return alltimes


#now loop through and find the neighbours for everything

#define column labels
columnlabels = ['Metadata_FrameNumber', 'Leading Cell ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']

print "Finding neighbours of leading edge red cells..."
np_neighbours_red = find_neighbours(reds, greens, red, green)
np_neighbours_red = pd.DataFrame(np_neighbours_red, columns=columnlabels)

print "Finding neighbours of leading edge green cells..."
np_neighbours_green = find_neighbours(greens, reds, green, red)    
np_neighbours_green = pd.DataFrame(np_neighbours_green, columns=columnlabels)

green_csv_name = 'allneighbours' + '_green_1.csv'
red_csv_name = 'allneighbours' + '_red_1.csv'

count = 1
while os.path.isfile(red_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    red_csv_name = 'leadingedge_red' + str(count) + '.csv'

count = 1
while os.path.isfile(green_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    green_csv_name = 'leadingedge_green' + str(count) + '.csv'

np_neighbours_red.to_csv(red_csv_name, index=False, header=True, sep=',')
np_neighbours_green.to_csv(green_csv_name, index=False, header=True, sep=',')

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print 'CSVs produced: ' + green_csv_name + ' and ' + red_csv_name
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'

