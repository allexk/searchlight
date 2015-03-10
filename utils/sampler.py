#!/usr/bin/python

import argparse
import csv
import cStringIO as StringIO
import os
import subprocess

def get_meta(iquery, array, what):
    """Retrieves metadata from a SciDb array.

    The metadata is presented as an array of dicts, one dict per
    object (attribute, dimension, etc.)

    Parameters:
        iquery -- path to the iquery program
        array  -- the name of the array to query
        what   -- what to query: dims, attrs
    """
    # get rid of any version ids
    ver_pos = array.rfind('@')
    if ver_pos != -1:
        array_nover = array[:ver_pos]
    else:
        array_nover = array

    # we want AFL and CSV output with a header
    args = [iquery, '-a', '-o', 'csv+']
    if what == 'dims':
        query = 'dimensions(%s)' % array_nover
    elif what == 'attrs':
        query = 'attributes(%s)' % array_nover
    else:
        raise ValueError, 'Unknown type of metadata: %s' % what
    args.extend(['-q', query])

    # the output (should be CSV with a header)
    out = StringIO.StringIO(subprocess.check_output(args,
        stderr=subprocess.STDOUT))
    csv_reader = csv.DictReader(out, quotechar="'")
    res = []
    for row in csv_reader:
        res.append(row)

    return res

def sample_attribute(iquery, array, chunks, attr, region):
    """ Samples the specified attribute from the array.

    The result is returned as a StringIO object.

    Parameters:
        iquery       -- path to the iquery program
        array        -- the name of the array to query
        chunks       -- a list of chunk sizes to use
        attr         -- the attribute to sample
        aid          -- the attribute's id
        total_chunks -- the total number of sample chunks
        region      -- particular region to sample
    """
    args = [iquery, '-a', '-o', 'csv+', '-q']
    array_spec = array
    if region:
        region_lbs = ', '.join([str(i) for i in region[0]])
        region_rbs = ', '.join([str(i) for i in region[1]])
        array_spec = 'between(%s, %s, %s)' % (array, region_lbs, region_rbs)

    chunks_spec = ', '.join([str(i) for i in chunks])
    regrid_cmd = 'regrid(%s, %s, ' % (array_spec, chunks_spec)
    regrid_cmd += 'min(%s) as min, max(%s) as max, sum(%s) as sum, '\
        'count(%s) as count)' % ((attr, ) * 4)
    print 'Sampling with the command: %s' % regrid_cmd

    # the full command
    args.append(regrid_cmd)
    out = StringIO.StringIO(subprocess.check_output(args,
        stderr=subprocess.STDOUT))

    return out


# command line parser
parser = argparse.ArgumentParser(description='''Samples the specified SciDb
    array, dumps the sample into a CSV file and creates a script for loading''')

# command line parameters
parser.add_argument('--chunks', metavar='N N...', nargs='+', type=int,
    help='Sample chunk sizes')
parser.add_argument('--bindir', metavar='path', default=os.getcwd(),
    help='Path to the SciDb binaries')
parser.add_argument('--outcsv', metavar='filepath', default=None,
    help='File to write the sample CSV')
parser.add_argument('--outsh', metavar='filepath', default=None,
    help='File to write the shell script for loading the sample from CSV')
parser.add_argument('--region', metavar='N', default=None, nargs='+', type=int,
    help='Specific region to sample: [lbs1, lbs2, ... rbs1, rbs2, ...], ...'
        '(inclusive, must be aligned with the chunks)')
parser.add_argument('array', help='Array to sample')
parser.add_argument('attr', help='Attribute to sample')

# parse
opts = parser.parse_args()

# full iquery path
iquery_path = os.path.join(opts.bindir, 'iquery')

# get dimensions and attributes
dims = get_meta(iquery_path, opts.array, 'dims')
attrs = get_meta(iquery_path, opts.array, 'attrs')

# one chunk size per dimension
if len(dims) != len(opts.chunks):
    raise ValueError, "The number of chunk sizes should be %d" % len(dims)

# check the area
regions = None
if opts.region:
    regions = []
    if len(opts.region) % (2 * len(dims)) != 0:
        raise ValueError, 'Sample regions must have %d dimensions' % len(dims)
    for i in range(0, len(opts.region), 2 * len(dims)):
        region_lbs = opts.region[i:i + len(dims)]
        region_rbs = opts.region[i + len(dims):i + 2 * len(dims)]
        for i in range(len(region_lbs)):
            chunk_i = opts.chunks[i]
            len_i = region_rbs[i] - region_lbs[i] + 1
            if region_lbs[i] % chunk_i != 0 or len_i % chunk_i != 0:
                raise ValueError, 'Area bounds must be aligned with the chunk!'
        regions.append((region_lbs, region_rbs))
    print 'Detected the following sample regions: %s' % str(regions)

# determine the name of the sample array and the total chunks number
sample_array_name = opts.array + '_' + opts.attr + '_'
array_origin = []
sample_array_end = []
for (i, chunk) in enumerate(opts.chunks):
    if i != 0:
        sample_array_name += 'x'
    sample_array_name += str(chunk)

    array_origin.append(int(dims[i]['start']))
    sample_array_end.append((int(dims[i]['length']) - 1) / chunk)

# sample each attribute into the file
csv_file_name = opts.outcsv
if not csv_file_name:
    csv_file_name = '%s.csv' % sample_array_name
# let's get the full name for the load script later
csv_file_name = os.path.abspath(csv_file_name)
csv_file = open(csv_file_name, 'w')

samples = []
header = None
if regions:
    for r in regions:
        sample = sample_attribute(iquery_path, opts.array, opts.chunks,
                                 opts.attr, r)
        # header is the same and we remove it from every StringIO "file"
        header = sample.readline().strip()
        samples.append(sample)
else:
    sample = sample_attribute(iquery_path, opts.array, opts.chunks,
                             opts.attr, None)
    header = sample.readline().strip()
    samples.append(sample)

# write to the file
csv_file.write(header)
for s in samples:
    for l in s:
        # We need to convert coordinates to the synopsis one here.
        # The way regrid() works is by setting the first coordinate
        # equivalent to the original coordinate and then incrementing them
        # when going to the next chunk.
        l = l.strip().split(',')
        for i in range(len(dims)):
            l[i] = str(int(l[i]) - array_origin[i])
        csv_file.write('\n%s' % ','.join(l))
csv_file.close()

# create the script
script_file_name = opts.outsh
if not script_file_name:
    script_file_name = 'load_%s.sh' % sample_array_name
script_file = open(script_file_name, 'w')

raw_sample_array_name = 'raw_' + sample_array_name
script_file.write('#!/bin/bash\n\n')

script_file.write('# Exit on error\nset -e\n\n')

script_file.write('# Set the proper PATH\n')
script_file.write("PATH=%s:$PATH\n\n" % opts.bindir)

script_file.write('# Remove old and create a new raw array\n')
script_file.write("iquery --ignore-errors -a -q 'remove(%s)'\n" % \
    raw_sample_array_name)
script_file.write("iquery -q 'create array %s <" % raw_sample_array_name)
for d in dims:
    script_file.write('%s: int64, ' % d['name'])
script_file.write("min: double, max: double, sum: double, count: uint64>\
[i=0:*,10000,0]'\n\n")

sdb_file_name = os.path.join(os.path.dirname(csv_file_name),
        '%s.sdb' % sample_array_name)
script_file.write('# Converting to scidb format\n')
script_file.write("csv2scidb -i %s -o %s -c 10000 -f 0 -p NNNNNN -s 1\n\n" % \
    (csv_file_name, sdb_file_name))

script_file.write('# Loading the CSV into the raw array\n')
script_file.write("iquery -n -q 'load %s from '\\''%s'\\'\n\n" % \
    (raw_sample_array_name, sdb_file_name))

script_file.write('# Creating the sample array\n')
script_file.write("iquery -q 'create array %s <min: double, max: double, \
sum: double, count: uint64>[" % sample_array_name)
for (i, d) in enumerate(dims):
    if i > 0:
        script_file.write(', ')
    script_file.write('%s=0:%d,1000,0' % (d['name'], sample_array_end[i]))
script_file.write("]'\n\n")

script_file.write('# Converting raw to sample\n')
script_file.write("iquery -n -a -q 'store(redimension(%s, %s), %s)'\n" % \
    (raw_sample_array_name, sample_array_name, sample_array_name))

script_file.close()
