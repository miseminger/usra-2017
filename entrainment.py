#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 31 May 2017

Processes a simulation xml file.
Arguments (1):
    full path to file
Outputs:
    creates graph of target parameter of all cells vs. time and saves in current directory

cross correlation: graph dot product vs. data shift and find maxima
"""

from lxml import etree as et
import numpy as np
import sys
import matplotlib.pyplot as plt
import math

if len(sys.argv) == 1:
    print '\nError: no file path specified\n'
    exit()

tree = et.parse(sys.argv[1])
root = tree.getroot()
timeframe_elements = root.findall(".//time")

# Building data array
data = 0  # rows = timeframes, columns = cells
timeframe_num = len(timeframe_elements)
cell_num = len(timeframe_elements[0].findall(".//cell"))

for timeframe in timeframe_elements:
    cell_elements = timeframe.findall(".//cell")
    timeframe_data = np.empty(cell_num)

    for cell in cell_elements:
        timeframe_data[cell_elements.index(cell)] = cell.get('gtpase')
    if timeframe_elements.index(timeframe) == 0:
        data = [timeframe_data]
    else:
        data = np.concatenate((data, [timeframe_data]))

# Plot

plt.plot(data)
plt.xlabel('time')
plt.ylabel('gtpase')
plt.title('gtpase vs. time')
plt.grid(True)
plt.savefig('gtpase_plot.png')

# Cross Correlation with first cell (index 0, cell_id="1")

data = np.transpose(data)  # rows = cells, columns = timeframes
ref_index = 0  # index of reference cell

# opt_shift
# Arguments:
#   ref_data - reference data set
#   shift_data - data to be compared to reference
# Output:
#   array of cross correlated data (dot product after shifting)

def opt_shift(ref_data, shift_data):
    shift_product = np.empty(len(shift_data))

    shift = - int(math.ceil(len(shift_data) / 2.0))
    while shift < math.floor(len(shift_data) / 2.0):
        shift_product[shift] = np.inner(ref_data, np.roll(shift_data, shift))

        shift += 1

    return shift_product

shift_products = np.empty([cell_num, timeframe_num])
cross = np.empty(cell_num)

cell_index = 0
while cell_index < cell_num:
    shift_products[cell_index] = opt_shift(data[ref_index], data[cell_index])
    cross[cell_index] = int(np.where(shift_products[cell_index] == max(shift_products[cell_index]))[0][0])
    cell_index += 1

def test():
    t = np.arange(0, 2, 0.004)
    a = np.empty(500)
    b = np.empty(500)

    i = 0
    while i < 500:
        a[i] = math.sin(2 * math.pi * t[i]) + np.random.normal(0, 0.5)
        b[i] = math.sin(2 * math.pi * (t[i] - 1.2)) + np.random.normal(0, 0.5)

        i += 1

    plt.plot(np.arange(-1, 1, 0.004), a, np.arange(-1, 1, 0.004), b)
    plt.show()

    plt.plot(np.arange(-1, 1, 0.004), opt_shift(a, b))
    plt.show()

plt.plot(np.transpose(shift_products))
plt.savefig('gtpase_shift.png')
