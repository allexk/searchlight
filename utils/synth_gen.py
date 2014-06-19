# -*- coding: utf-8 -*-
# Creates synthetic data according to the parameters and dumps it into a CSV

import numpy as np
import argparse
import random
import cStringIO as StringIO
import csv
import math

def row_major_iter(lbs, rbs, steps):
    """Generates row-major traversal of the specified array.

    The array is traversed from lbs to rbs (inclusive) with
    strides specified in steps.

    Parameters:
        lbs -- the leftmost corner
        rbs -- the rightmost corner
        steps -- the stride sizes
    """
    point = list(lbs)
    dims = len(lbs)
    have_more_points = True
    while have_more_points:
        yield point
        # move to the next one
        i = dims - 1
        while True:
            point[i] += steps[i]
            if point[i] > rbs[i]:
                if i == 0:
                    have_more_points = False
                    break
                else:
                    point[i] = lbs[i]
                    i -= 1
            else:
                break

def get_rbs_for_lbs(lbs, lens, rbounds):
    """Computes the right boundary of a region for the specified left boundary.

    The boundary is computed according to the specified steps. The
    boundary cannot go beyond the specified on in the parameters. Otherwise,
    it is truncated.

    Parameters:
        lbs -- the left boundary of the region
        lens -- the length of the region
        rbounds -- the maximum right boundary
    """
    dims = len(lbs)
    rbs = []
    for i in range(dims):
        rbs_i = lbs[i] + lens[i] - 1
        if rbs_i > rbounds[i]:
            rbs_i = rbounds[i]
        rbs.append(rbs_i)
    return rbs


# parse args
parser = argparse.ArgumentParser(description='Creates synthetic data'
                        'according to the parameters and dumps it into a CSV')
parser.add_argument('--lbs', metavar='N', nargs='+', type=int, required=True,
                    help='Specifies the leftmost corner')
parser.add_argument('--rbs', metavar='N', nargs='+', type=int, required=True,
                    help='Specifies the rightmost corner')
parser.add_argument('--sgrid', metavar='N', nargs='+', type=int, required=True,
                    help='Specifies the small grid of cells')
parser.add_argument('--count', metavar='N', default=0, type=long,
                    help='Specifies the number of points')
parser.add_argument('--density', metavar='N.N', default=1.0, type=float,
                    help = 'Specifies the density of the array')
parser.add_argument('--lgrid', metavar='N', default=[], nargs='*', type=int,
                    help = 'Specifies the large (cluster) grid')
parser.add_argument('--clusts', metavar='N', default=[0], nargs='*',
                    help = 'Specifies the sub-areas to populate')
parser.add_argument('--mean', metavar='N', default=[50, 201, 25],
                    nargs=3, type=int,
                    help='Means for the normal generator: [min, max) + step')
parser.add_argument('--dev', metavar='N', default=5.0, type=float,
                    help='Standard deviation for the generator')
parser.add_argument('file', help='Filename for the output')
opts = parser.parse_args()

# make proper means
means = list(range(opts.mean[0], opts.mean[1], opts.mean[2]))
print 'Using the following means: %s' % str(means)
print 'Using the following deviation: %.3f' % opts.dev

# make some conversions
lbs = opts.lbs
rbs = opts.rbs
sgrid = opts.sgrid
dims = len(lbs)
assert len(lbs) == len(rbs) == len(sgrid) != 0, \
    'lbs, rbs and step specifications must match!'

# the total number of elements
total_elems = 1
for i in range(dims):
    total_elems *= rbs[i] - lbs[1] + 1
density = opts.density
if opts.count != 0:
    density = float(opts.count) / total_elems

# large grid (specifies the large clusters)
lgrid = opts.lgrid
if len(lgrid) == 0:
    for i in range(dims):
        lgrid.append(rbs[i] - lbs[i] + 1)
    print 'The large grid is not specified, default=%s' % str(lgrid)
else:
    assert len(lbs) == len(lgrid), \
        'The large grid must have the same dimensionality!'

clusters = []
cluster_generator = row_major_iter(lbs, rbs, lgrid)
for clust_lbs in cluster_generator:
    clust_rbs = get_rbs_for_lbs(clust_lbs, lgrid, rbs)
    clusters.append((list(clust_lbs), clust_rbs))
# print
print 'Working with the following clusters: %s' % str(clusters)

csv_file = open(opts.file, 'wb')
for i in range(dims):
    csv_file.write('x%d,' % i)
csv_file.write('val\n')
csv_writer = csv.writer(csv_file)

for clust_num in opts.clusts:
    cluster = clusters[clust_num]
    clust_lbs = cluster[0]
    clust_rbs = cluster[1]

    print 'Generating tuples for area %s with density %.3f' % (str(cluster), density)
    cell_gen = row_major_iter(clust_lbs, clust_rbs, sgrid)
    for cell_lbs in cell_gen:
        cell_rbs = get_rbs_for_lbs(cell_lbs, sgrid, clust_rbs)
        # iterate over the cell and generate
        cell_tuples = []
        cell_mean = means[np.random.randint(len(means))]
        elem_gen = row_major_iter(cell_lbs, cell_rbs, [1] * dims)
        for elem in elem_gen:
            if np.random.random_sample() <= density:
                val = round(np.random.normal(cell_mean, opts.dev), 3)
                cell_tuples.append(elem + [val])
        # write to the file
        csv_writer.writerows(cell_tuples)

csv_file.close()

print "Finished!"
