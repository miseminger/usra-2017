#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 31 May 2017

Extracts information from a metadata xml file.
Arguments (1):
    full path to file
Outputs:
    dimensions (X,Y) in pixels and microns
    time (T) in # frames and minutes
    channel number corresponding to each LUT name
"""

from lxml import etree as et
import numpy as np
import re
import sys

if len(sys.argv) == 1:
    print '\nError: no file path specified\n'
    exit()

tree = et.parse(sys.argv[1])
root = tree.getroot()

dimension_elements = root.findall(".//DimensionDescription")
dimension_values = np.empty([3, 2])
print '\nDimensions:'

for dim in dimension_elements:
    if dim.get("DimID") == "T" or dim.get("DimID") == "4":
        print '{:3} {:>10} {:10}'.format('T', dim.get("NumberOfElements"), 'frames'),
        if dim.get("Unit") == "s":
            time = '%.3e' % (float(dim.get("Length")))
            unit = "seconds"
            frames_per_second = float(dim.get("NumberOfElements")) / float(time)
        else:
            time_array = np.delete(np.array(re.split('[hms]', dimension_elements[2].get("Length"))), 3).astype(float)
            time = 60 * time_array[0] + time_array[1] + time_array[2] / 60
            unit = "minutes"
            frames_per_second = float(dim.get("NumberOfElements")) / (time * 60)
        print '{:>10} {:10}'.format(time, unit)
    else:
        if dim.get("Unit") == "m":
            length = float(dim.get("Length")) * 10 ** 6
        else:
            length = dim.get("Length")
        if dim.get("DimID") == "X" or dim.get("DimID") == "1":
            dim_name = "X"
            microns_per_pixel = float(length) / float(dim.get("NumberOfElements"))
        else:
            dim_name = "Y"
        print '{:3} {:>10} {:10} {:>10} {:10}'.format(dim_name, dim.get("NumberOfElements"), 'pixels', length,
                                                      'microns')

wfci_elements = root.findall(".//WideFieldChannelInfo")
wfci_values = np.empty(3)
print '{:7} {:15}'.format('\nLUT', 'Channel Number')

for wcfi in wfci_elements:
    print '{:6} {:15}'.format(wcfi.get("LUT"), wcfi.get("Channel")[-2:])

print '{:20} {:>5.4} {:20} {:>5.4}'.format('\nMicrons per pixel:', microns_per_pixel, '\nFrames per second:',
                                           frames_per_second)
print ''
