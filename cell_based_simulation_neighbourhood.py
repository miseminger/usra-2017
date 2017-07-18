#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 08 June 2017

Processes a Compucell3D simulation xml file and constructs a graph of all the cells at a given timeframe by location and
indicates neighbouring cells with a connection.
Arguments:
    -f inputfile - full path to XML file
    -t timeframe - time at which cell neighbour graph is generated
Output: saves image of graph in the script directory.
"""

import math

from lxml import etree as et
from optparse import OptionParser
from graphviz import Graph

parser = OptionParser()
parser.add_option('-f', '--file', type='string', dest='inputfile', help='path to XML file')
parser.add_option('-t', '--time', type='int', dest='time', help='timeframe to graph')
parser.add_option('-d', '--printdist', action='store_true', dest='printdist', help='print distance on edges')

(options, args) = parser.parse_args()
inputfile = options.inputfile
time = options.time
printdist = 0
if options.printdist:
    printdist = options.printdist

tree = et.parse(inputfile)
root = tree.getroot()

timeframe = root[time]

redf = '#FF9090'
greenf = '#75CF90'
redl = '#FF3535'
greenl = '#019A2F'

g = Graph(name='cell_neighbours_t_' + str(time), format='png', engine='neato')

g.attr('node', pin='true', shape='circle', style='filled', width='0.5', fixedsize='true', penwidth='0',
       fontname='courier')
g.attr('edge', fontsize='10.0', penwidth='4', fontname='courier')

cell_index = 0
while cell_index < len(timeframe):
    cell_type = timeframe[cell_index].get('type')
    color = greenf
    if int(cell_type) == 2:
        color = redf
    x = float(timeframe[cell_index].get('x')) / 20
    y = float(timeframe[cell_index].get('y')) / 20
    position = str(x) + ',' + str(y)
    neighbors = timeframe[cell_index].get('neighbors')
    if neighbors != '':
        neighbors = map(int, neighbors.split(' '))

    g.node(str(cell_index), pos=position, fillcolor=color)

    for neighbor in neighbors:
        if cell_index < neighbor:
            if cell_type == timeframe[neighbor].get('type'):
                lcolor = color
            else:
                lcolor = 'black'
            nx = float(timeframe[neighbor].get('x'))
            ny = float(timeframe[neighbor].get('y'))
            distance = ''
            if printdist:
                distance = str(round(math.sqrt((x - nx)**2 + (y - ny)**2), 2))

            g.edge(str(cell_index), str(neighbor), color=lcolor, label=distance)

    cell_index += 1

g.render()
