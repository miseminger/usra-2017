#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 6 June 2017

Processes a Compucell3D simulation xml file and graphs one or more parameters.
Arguments:
    -f inputfile - full path to XML file
    -p parameters - comma separated parameters to graph
    -i image - method of graphing multiple parameters
        0 (default) - plot all parameters on the same axis, one output
        1 - plot each parameter on its own subplot of a figure, one output of the entire figure
        2 - plot each paramater on its own figure, one output for each parameter
Output: graphs max, min, and average of a parameter and saves image(s) of graph in the script directory.
"""

import math
import numpy as np
import matplotlib.pyplot as plt

from lxml import etree as et
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--file', type='string', dest='inputfile', help='path to XML file')
parser.add_option('-p', '--parameters', type='string', dest='parameters', help='comma separated parameters to graph')
parser.add_option('-i', '--image', type='int', dest='image', help='0 (default) - plot on same axis, one output; ' +
        '1 - plot separately, one output with all plots; 2 - plot separately, one output for each plot')

(options, args) = parser.parse_args()
inputfile = options.inputfile
parameters = options.parameters.split(',')
image = 0
if options.image:
    image = options.image

tree = et.parse(inputfile)
root = tree.getroot()

data = np.empty([len(parameters), len(root), 3])

parameter = 0
while parameter < len(parameters):

    time = 0
    while time < len(root):
        parameter_values = np.empty(len(root[time]))

        cell = 0
        while cell < len(root[time]):
            parameter_values[cell] = root[time][cell].get(parameters[parameter])
            cell += 1

        data[parameter, time, 0], data[parameter, time, 1], data[parameter, time, 2] = \
            np.amax(parameter_values), np.amin(parameter_values), np.average(parameter_values)

        time += 1
    parameter += 1

if len(parameters) == 1:
    plt.plot(data[0])
    plt.legend(['max', 'min', 'average'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=9)
    plt.ylabel(parameters[0])
    plt.xlabel('time')
    plt.show()
else:
    if image == 0:
        plt.subplot2grid((1, 4), (0, 0), colspan=3)

        parameter = 0
        while parameter < len(parameters):
            plt.plot(np.transpose(data[parameter])[0], label=parameters[parameter] + ' max')
            plt.plot(np.transpose(data[parameter])[1], label=parameters[parameter] + ' min')
            plt.plot(np.transpose(data[parameter])[2], label=parameters[parameter] + ' average')

            parameter += 1
        plt.xlabel('time')
        plt.ylabel(', '.join(parameters))
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=9)
        plt.savefig('_'.join(parameters) + '.png')
    if image == 1:
        side = int(math.ceil(math.sqrt(len(parameters))))
        maxrows = int(math.ceil(float(len(parameters)) / side))

        plt.figure(1)

        parameter = 0
        while parameter < len(parameters):
            plt.subplot(maxrows, side, parameter + 1)

            plt.plot(data[parameter])
            plt.xlabel('time')
            plt.ylabel(parameters[parameter])

            parameter += 1
        plt.tight_layout()
        plt.savefig('_'.join(parameters) + '.png')
    if image == 2:
        parameter = 0
        while parameter < len(parameters):
            plt.subplot2grid((1, 4), (0, 0), colspan=3)

            plt.plot(np.transpose(data[parameter])[0], label=parameters[parameter] + ' max')
            plt.plot(np.transpose(data[parameter])[1], label=parameters[parameter] + ' min')
            plt.plot(np.transpose(data[parameter])[2], label=parameters[parameter] + ' average')

            plt.xlabel('time')
            plt.ylabel(parameters[parameter])
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=9)

            plt.show()

            parameter += 1
