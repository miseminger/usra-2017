#!/usr/bin/env python

"""
Author: Cindy Tan
Last Modified Date: 3 August 2017

OPTICS ordering implementation in python
Method to extract clusters from reachability plot

Usage: Optics(dataset).cluster()
   - dataset is an array of Point instances:
       - initialize with Point(*coordinates)
   - returns Clustering instance with attributes:
       - .ordered: OPTICS ordering of points
       - .labels: clustering labels of ordered
       - .centroids: cluster centroids
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt

class Point:

    def __init__(self, *coordinates):
        self.pos = coordinates
        self.dim = len(coordinates)

    def __repr__(self):
        return str(self.pos)

class Clustering:

    def __init__(self, ordered_points, labels, centroids):
        self.ordered = ordered_points
        self.labels = labels
        self.centroids = centroids

class Optics:

    def __init__(self, points, eps, min_pts):
        self.unprocessed = points
        self.eps = eps
        self.min_pts = min_pts

    def _setup(self):
        for p in self.unprocessed:
            p.reach_dist = None
        self.ordered = []

    def _distance(self, p, q):
        sum = 0
        i = 0
        while i < p.dim:
            sum += (p.pos[i] - q.pos[i])**2
            i += 1
        return math.sqrt(sum)

    def _neighbours(self, p):
        neighbours = []
        for q in self.unprocessed:
            if q is not p and self._distance(p, q) <= self.eps:
                neighbours.append(q)
        return neighbours

    def _process(self, p):
        self.unprocessed.remove(p)
        self.ordered.append(p)
        if len(self._neighbours(p)) >= self.min_pts - 1:
            neighbour_dists = sorted([self._distance(p, n) for n in self._neighbours(p)])
            p.core_dist = neighbour_dists[self.min_pts - 2]
        else:
            p.core_dist = None

    def _update(self, p, queue):
        for n in [n for n in self._neighbours(p) if n in self.unprocessed]:
            reach_dist = max(p.core_dist, self._distance(p, n))
            if n.reach_dist is None:
                n.reach_dist = reach_dist
                queue.append(n)
            elif reach_dist < n.reach_dist:
                n.reach_dist = reach_dist

    def order(self):

        # Initialize all reachability distances to be None
        # Create 'unprocessed' list of points
        # Create empty list for output of ordered points
        self._setup()

        # For each unprocessed point p:
        #   - remove from 'unprocessed' list
        #   - add to ordered list
        #   - determine a core distance
        while self.unprocessed:
            p = self.unprocessed[0]
            self._process(p)

            # If p is a core point, create 'queue' to process its neighbours
            if p.core_dist is not None:
                queue = []

                # Compute reachability distances of all neighbours of p and add neighbours to 'queue'
                self._update(p, queue)

                # Until the queue is empty:
                while queue:

                    # Sort queue by reachability distance, remove first point
                    queue.sort(key=lambda n: n.reach_dist)
                    n = queue.pop(0)

                    # For neighbour n:
                    #   - remove from 'unprocessed' list
                    #   - add to ordered list
                    #   - determine a core distance
                    self._process(n)

                    # If n is a core point, compute reachability distances of all neighbours of n
                    if n.core_dist is not None:
                        self._update(n, queue)

        return self.ordered

    def cluster(self):
        
        # Create directory and save figures
        directory_i = 1
        directory = 'OPTICS_images/F_' + str(directory_i)
        while os.path.exists(directory):
            directory_i += 1
            directory = 'OPTICS_images/F_' + str(directory_i)
        os.makedirs(directory)
        
        def savefig(name):
            count = 1
            while os.path.isfile(directory + '/' + name + '_' + str(count) + '.png'):
                count += 1
            plt.savefig(directory + '/' + name + '_' + str(count) + '.png')
        
        # Zero-centred standard deviation
        def std(data):
            return math.sqrt(np.mean([d**2 for d in data]))
            
            
        ## ARRAY DEFINITIONS
        
        # Reachability distances
        ordered = self.order()
        reach_dists = [p.reach_dist for p in ordered]
        reach_dists[0] = reach_dists[1]

        # Derivative
        derivative = []
        blur_radius = 1  # Radius around point to take discrete derivative
        for i in range(len(reach_dists)):
            
            # Construct radius using j, if full radius cannot be constructed, use an asymmetrically reduced radius as necessary
            j = i - blur_radius
            count = 0
            minfound = False
            while j <= i + blur_radius:
                if j >= 0:
                    try:
                        maxval = reach_dists[j]
                        count += 1
                        if not minfound:
                            minval = reach_dists[j]
                            minfound = True
                    except IndexError:
                        pass
                j += 1
            derivative.append((maxval - minval) / count)
        
        # Plot all arrays
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(reach_dists)
        axarr[1].plot(derivative)
        axarr[1].axhline(std(derivative), c='black')  # Plot threshold values (standard deviation)
        axarr[1].axhline(-std(derivative), c='black')
        savefig('data_arrays')
        
        
        ## CREATE PARTITIONS
            
        # Find dividers by finding the approximate noise between clusters
        dividers = []         # Array of dividers
        divider_start = None  # Start of the dividing noise region
        divider_end = None    # End of the dividing noise region
        noise_estimate = 0    # Total count of points in all noise regions
        
        # Given the start and end of a noise region, finds the exact dividing point at the centre
        def divider(start, end):
            divider = float(start + end) / 2
            
            # If centre of noise region is a single point, define it as the divider
            if int(divider) == divider:
                divider = int(divider)
                
            # Otherwise, of two indices around the centre, whichever has smaller derivative magnitude is the divider
            else: 
                if abs(derivative[int(math.ceil(divider))]) <= abs(derivative[int(math.floor(divider))]):
                    divider = int(math.ceil(divider))
                else:
                    divider = int(math.floor(divider))
            return divider
        
        # Iterate through derivative array
        for i in range(len(derivative)):
            
            # If derivative is larger in magnitude than a standard deviation, do something:
            # If derivative is positive, start a noise region
            if derivative[i] >= std(derivative):
                
                # If a noise region is already being worked on and it has an endpoint, finish that noise region and create a divider
                # If it has no endpoint, do nothing
                if divider_end is not None:
                    dividers.append(divider(divider_start, divider_end))
                    noise_estimate += divider_end - divider_start
                    divider_end = None
                    divider_start = None
                
                # If no noise regions are currently being worked on, start a new noise region
                if divider_start is None:
                    divider_start = i
            
            # If derivative is negative and there is a noise region being worked on, give it an endpoint
            # If the noise region already has an endpoint, replace the old endpoint with the current point
            if derivative[i] <= -std(derivative) and divider_start is not None:
                divider_end = i - 1
            
            # If the end of the array is reached
            if i == len(derivative) - 1:
                
                # If the current noise region has an endpoint, finish the noise region and create a divider
                if divider_end is not None:
                    dividers.append(divider(divider_start, divider_end))
                    noise_estimate += divider_end - divider_start
                
                # If the current noise region has no endpoint, count the region as noise but do not create a divider
                elif divider_start is not None:
                    noise_estimate += len(derivative) - divider_start
        
        partitions = np.split(np.array(reach_dists), dividers)  # Create partitions of reach_dists array using dividers
        noise_estimate = float(noise_estimate) / len(ordered)   # Quantify noise as a percentage
        
        
        ## COMPLETE CLUSTERS
        
        # Get clusters
        clusters = []                 # Array of clusters
        noise = [p for p in ordered]  # Array of noise (copy of all points, non-noise points will be removed later)
        labels = [-1] * len(ordered)  # Array of cluster labels for each point (default to -1, noise cluster)
        centroids = []                # Array of cluster centroids
        points_processed = 0          # Number of points processed
        f, axes = plt.subplots(1, len(partitions), sharey=True)  # Start a figure with subplots for partitions
        
        # If noise is estimated to be less than 5%, accept partitions as complete clusters with no noise
        if noise_estimate <= 0.05:
            for i in range(len(partitions)):
                partition = partitions[i]
                
                cluster_points = ordered[points_processed:points_processed + len(partition)]
                points_processed += len(partition)
                for p in cluster_points:
                    noise.remove(p)
                clusters.append(Cluster(*cluster_points))
                
                cluster_id = len(clusters)
                labels[points_processed:points_processed + len(partition)] = [cluster_id] * len(cluster_points)
                centroid.append(Point(*[s / len(cluster_points) for s in sums]))
                
                axes[i].plot(partition)
        
        # Otherwise, if noise estimated to be significant, find cluster endpoints
        else:  
            
            # Given data and a vertical threshold, define cluster in partition as all points below threshold
            def process(reach_data, threshold):
                start, end = None, None
                for j in range(len(reach_data)):
                    if reach_data[j] >= threshold and start is not None and end is None:
                        end = min(j + blur_radius, len(reach_data))
                    if reach_data[j] <= threshold and start is None:
                        start = max(j - blur_radius, 0)
                if end is None:
                    end = len(reach_data)
                axes[i].axhline(threshold)
                return start, end
            
            # Iterate through partitions
            for i in range(len(partitions)):
                partition = partitions[i]
                axes[i].plot(partition)
                
                # Convergent recursion to process points:
                #   - start by processing with threshold 2 * minimum value of partition
                #   - for the next iteration, find the value of the centre point of the past iteration and use 2 times that value as the threshold
                #   - iterate until there is no change in cluster endpoints
                old_cluster_endpoints = (0, len(partition))
                cluster_endpoints = process(partition, 2 * min(partition))
                count = 0
                while cluster_endpoints != old_cluster_endpoints and count < 10:
                    old_cluster_endpoints = cluster_endpoints
                    cluster_endpoints = process(partition, 2 * partition[int(math.ceil(np.mean(cluster_endpoints)))])
                    count += 1
                
                # Convert from partition indexing (starts from 0 for every partition) to all data indexing
                cluster_endpoints = [cluster_endpoints[0] + points_processed, cluster_endpoints[1] + points_processed]
                points_processed += len(partition)
                
                # Get points, create Cluster instance, and remove from noise
                cluster_points = ordered[cluster_endpoints[0]:cluster_endpoints[1]]
                for p in cluster_points:
                    sums = np.zeros(p.dim)
                    for c in range(p.dim):
                        sums[c] += p.pos[c]
                    noise.remove(p)
                
                cluster_id = len(clusters)
                labels[cluster_endpoints[0]:cluster_endpoints[1]] = [cluster_id] * len(cluster_points)
                centroids.append(Point(*[s / len(cluster_points) for s in sums]))
                
                clusters.append(cluster_points)
        
        clustering = Clustering(ordered, labels, centroids)
                
        savefig('partitions')
        
        # Plot clusters
        # colourmap   red        blue       yellow     pink       green      orange     purple     sky blue
        colormap = ['#FF0033', '#0033FF', '#FFFF33', '#FF0099', '#33CC00', '#FF9900', '#9900FF', '#03A9F4',
                    # blue-grey  brown
                    '#34495E', '#9B4D00']
        
        plt.figure()
        for i in range(len(clusters)):
            plt.scatter([p.pos[0] for p in clusters[i]], [p.pos[1] for p in clusters[i]], c=colormap[i])
        plt.scatter([p.pos[0] for p in noise], [p.pos[1] for p in noise], c='w')
        savefig('clusters')
        
        return clustering
