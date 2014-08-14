#!/usr/bin/python

import argparse
import csv
import cStringIO as StringIO
import os
import subprocess
import math

def GetMeta(iquery, array, what):
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

def SampleAttribute(iquery, array, chunks, attr, aid, total_chunks, region):
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
    
    regrid_cmd += 'min(%s) as min, max(%s) as max, sum(%s) as sum,'\
        'count(%s) as count)' % ((attr, ) * 4)

    # reshape ignores attribute specs, but...
    reshape_cmd = 'reshape(' + regrid_cmd + ', <a: double>'

    # we need to specify proper dimensions
    reshape_cmd += '[chunk=0:%d,%d,0, attr=%d:%d,1,0])' % \
        (total_chunks - 1, total_chunks, aid, aid)
    print 'Sampling with the command: %s' % reshape_cmd

    # the full command
    args.append(reshape_cmd)
    out = StringIO.StringIO(subprocess.check_output(args,
        stderr=subprocess.STDOUT))

    return out


# command line parser
parser = argparse.ArgumentParser(description='''Samples the specified SciDb
    array, dumps the sample into a CSV file and creates a script for loading''')

# command line parameters
parser.add_argument('--attrs', metavar='att1 att2...', nargs='+',
    help='Attribute to sample')
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

# parse
opts = parser.parse_args()

# full iquery path
iquery_path = os.path.join(opts.bindir, 'iquery')

# get dimensions and attributes
dims = GetMeta(iquery_path, opts.array, 'dims')
attrs = GetMeta(iquery_path, opts.array, 'attrs')

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

attr_ids = {} # map of attribute names to ids
for oa in opts.attrs:
    found = False
    aid = -1
    for a in attrs:
        if a['name'] == oa:
            found = True
            aid = int(a['No'])
            break
    if not found:
        raise ValueError, "The array does not contain attribute %s" % oa
    attr_ids[oa] = aid

# determine the name of the sample array and the total chunks number
sample_array_name = opts.array + '_'
total_chunks = 1
for (i, chunk) in enumerate(opts.chunks):
    if i != 0:
        sample_array_name += 'x'
    sample_array_name += str(chunk)

    dim_len = float(dims[i]['length'])
    total_chunks *= math.ceil(dim_len / chunk)
total_chunks = int(total_chunks)

# sample each attribute into the file
csv_file_name = opts.outcsv
if not csv_file_name:
    csv_file_name = '%s.csv' % sample_array_name
# let's get the full name for the load script later
csv_file_name = os.path.abspath(csv_file_name)
csv_file = open(csv_file_name, 'w')

csv_header_needed = True
for (i, attr) in enumerate(opts.attrs):
    samples = []
    header = None
    if regions:
        for r in regions:
            sample = SampleAttribute(iquery_path, opts.array, opts.chunks,
                                     attr, attr_ids[attr], total_chunks, r)
            # header is the same and we remove it from every StringIO "file"
            header = sample.readline().strip()
            samples.append(sample)
    else:
        sample = SampleAttribute(iquery_path, opts.array, opts.chunks,
                                 attr, attr_ids[attr], total_chunks, None)
        header = sample.readline().strip()
        samples.append(sample)
    if csv_header_needed:
        csv_file.write(header)
        csv_header_needed = False
    for s in samples:
        for l in s:
            l = l.strip()
            csv_file.write('\n%s' % l)
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
script_file.write("iquery -q 'create array %s <chunk: int64, attr: int64, \
min: double, max: double, sum: double, count: uint64>[i=0:*,10000,0]'\n\n" % \
    raw_sample_array_name)

sdb_file_name = os.path.join(os.path.dirname(csv_file_name),
        '%s.sdb' % sample_array_name)
script_file.write('# Converting to scidb format\n')
script_file.write("csv2scidb -i %s -o %s -c 10000 -f 0 -p NNNNNN -s 1\n\n" % \
    (csv_file_name, sdb_file_name))

script_file.write('# Loading the CSV into the raw array\n')
script_file.write("iquery -n -q 'load %s from '\\''%s'\\'\n\n" % \
    (raw_sample_array_name, sdb_file_name))

script_file.write('# Creating the sample array\n')
script_file.write("iquery -q 'create array %s <min: double, \
max: double, sum: double, count: uint64>[chunk=0:%d,%d,0,\
attr=0:%d,%d,0]'\n\n" % (sample_array_name, total_chunks - 1,
        total_chunks, len(attrs) - 1, len(attrs)))

script_file.write('# Converting raw to sample\n')
script_file.write("iquery -n -a -q 'store(redimension(%s, %s), %s)'\n" % \
    (raw_sample_array_name, sample_array_name, sample_array_name))

script_file.close()
