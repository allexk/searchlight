#!/usr/bin/python

import argparse
import csv
import cStringIO as StringIO
import os


def get_meta(scidb_connection, array, what):
    """Retrieves metadata from a SciDb array.

    The metadata is presented as an array of dicts, one dict per
    object (attribute, dimension, etc.)

    Parameters:
        scidb_connection -- SciDB connection
        array  -- the name of the array to query
        what   -- what to query: dims, attrs
    """
    # get rid of any version ids
    ver_pos = array.rfind('@')
    if ver_pos != -1:
        array_nover = array[:ver_pos]
    else:
        array_nover = array
    # define query
    if what == 'dims':
        query = 'dimensions(%s)' % array_nover
    elif what == 'attrs':
        query = 'attributes(%s)' % array_nover
    else:
        raise ValueError, 'Unknown type of metadata: %s' % what
    # run via the connection
    query_id = scidb_connection.query(query, True, 'csv+')
    res_file = StringIO.StringIO(scidb_connection.get_result(query_id))
    csv_reader = csv.DictReader(res_file, quotechar="'")
    res = []
    for row in csv_reader:
        res.append(row)
    return res

def sample_attribute(scidb_connection, array, chunks, attr, region):
    """ Samples the specified attribute from the array.

    The result is returned as a StringIO object.

    Parameters:
        scidb_connection -- SciDB connection
        array        -- the name of the array to query
        chunks       -- a list of chunk sizes to use
        attr         -- the attribute to sample
        aid          -- the attribute's id
        total_chunks -- the total number of sample chunks
        region      -- particular region to sample
    """
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
    # run via the connection
    query_id = scidb_connection.query(regrid_cmd, True, 'csv+')
    res_file = StringIO.StringIO(scidb_connection.get_result(query_id))
    return res_file


# command line parser
parser = argparse.ArgumentParser(description='''Samples the specified SciDb
    array, dumps the sample into a CSV file and creates a script for loading''')

# command line parameters
parser.add_argument('--method', metavar='http/iquery', default='iquery',
                    help='Method to connect.')
parser.add_argument('--host', metavar='hostname', default='localhost',
                    help='Host to connect to')
parser.add_argument('--port', metavar='port', type=int, default=1239,
                    help='Port to connect to')
parser.add_argument('--user', metavar='name', default='SciDBUser',
                    help='User to connect as')
parser.add_argument('--password', metavar='path',
                    help='Path to a file with the password')
parser.add_argument('--bindir', metavar='path', default=os.getcwd(),
                    help='Path to the SciDb binaries')
parser.add_argument('--chunks', metavar='N N...', nargs='+', type=int,
                    help='Sample chunk sizes')
parser.add_argument('--outfile', metavar='filepath', default=None,
    help='File to write the sample')
parser.add_argument('--outscript', metavar='filepath', default=None,
    help='File to write the script for binary sample load')
parser.add_argument('--region', metavar='N', default=None, nargs='+', type=int,
    help='Specific region to sample: [lbs1, lbs2, ... rbs1, rbs2, ...], ...'
        '(inclusive, must be aligned with the chunks)')
parser.add_argument('--binary', action='store_true',
                    help='Write in the SciDB binary format')
parser.add_argument('array', help='Array to sample')
parser.add_argument('attr', help='Attribute to sample')

# parse
opts = parser.parse_args()

# create a connection
if opts.method == 'http':
    from scidb_http import SciDBConnection
    if not opts.password:
        raise ValueError('--password is mandatory for the http method')
    with open(opts.password, 'r') as passw_file:
        password = passw_file.read().strip()
    connection = SciDBConnection(opts.host, opts.port, opts.user, password)
    connection.connect()
elif opts.method == 'iquery':
    from scidb_iquery import SciDBConnection
    connection = SciDBConnection(opts.bindir, opts.port)
else:
    raise ValueError('Unknown connection method: %s' % opts.method)

# get dimensions and attributes
dims = get_meta(connection, opts.array, 'dims')

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
        for j in range(len(region_lbs)):
            chunk_j = opts.chunks[j]
            len_j = region_rbs[j] - region_lbs[j] + 1
            if region_lbs[j] % chunk_j != 0 or len_j % chunk_j != 0:
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
    # use start as total left boundary; needed for proper origin
    low_boundary = int(dims[i]['start'])
    high_boundary = int(dims[i]['high'])
    array_origin.append(low_boundary)
    curr_length = high_boundary - low_boundary + 1
    sample_array_end.append((curr_length - 1) / chunk)

# sample each attribute into the file
out_file_name = opts.outfile
if not out_file_name:
    if opts.binary:
        out_file_name = '%s.dat' % sample_array_name
    else:
        out_file_name = '%s.csv' % sample_array_name
# let's get the full name for the load script later
out_file_name = os.path.abspath(out_file_name)
out_file = open(out_file_name, 'wb')

samples = []
header = None
if regions:
    for r in regions:
        sample = sample_attribute(connection, opts.array, opts.chunks,
                                  opts.attr, r)
        # header is the same and we remove it from every StringIO "file"
        header = sample.readline().strip()
        samples.append(sample)
else:
    sample = sample_attribute(connection, opts.array, opts.chunks,
                              opts.attr, None)
    header = sample.readline().strip()
    samples.append(sample)

# write to the file
if opts.binary:
    import mimic
    binary_writer = mimic.BinarySciDBWriter(out_file, ['int64', 'int64',
                                                       'double null',
                                                       'double null',
                                                       'double', 'uint64'])
else:
    out_file.write(header)
for s in samples:
    for row in s:
        # We need to convert coordinates to the synopsis one here.
        # The way regrid() works is by setting the first coordinate
        # equivalent to the original coordinate and then incrementing them
        # when going to the next chunk.
        row = row.strip().split(',')
        for i in range(len(dims)):
            row[i] = str(long(row[i]) - array_origin[i])
        if opts.binary:
            row[0] = long(row[0])
            row[1] = long(row[1])
            if row[2] != 'null':
                row[2] = float(row[2])
            else:
                row[2] = None
            if row[3] != 'null':
                row[3] = float(row[3])
            else:
                row[3] = None
            row[4] = float(row[4])
            row[5] = long(row[5])
            binary_writer.writerow(row)
        else:
            out_file.write('\n%s' % ','.join(row))
out_file.close()

# create the script
# sample array schema
sample_array_schema = "<min: double null, max: double null, sum: double," \
                      "count: uint64>["
for (i, d) in enumerate(dims):
    if i > 0:
        sample_array_schema += ', '
    sample_array_schema += '%s=0:*,1000,0' % (d['name'])
sample_array_schema += ']'
# raw array schema (for input)
raw_array_schema = '<'
for d in dims:
    raw_array_schema += '%s: int64, ' % d['name']
raw_array_schema += "min: double null, max: double null, sum: double," \
                    "count: uint64>[i=0:*,10000,0]"
# binary tuple format
binary_tuple_schema = '(' + ', '.join(['int64'] * len(dims))
binary_tuple_schema += ', double null, double null, double, uint64)'

# create the script
script_file_name = opts.outscript
if not script_file_name:
    if opts.method == 'http':
        script_file_name = 'load_%s.script' % sample_array_name
    else:
        script_file_name = 'load_%s.sh' % sample_array_name
with open(script_file_name, 'w') as script_file:
    if opts.method == 'http':
        script_file.write('# create array (not needed if exists)\n')
        script_file.write('query:create array %s%s\n' %
                          (sample_array_name, sample_array_schema))

        script_file.write("# upload the binary to SciDB\n")
        script_file.write("upload:%s\n" % out_file_name)

        script_file.write("# load the sample (update if exists)\n")
        script_file.write("query:insert(redimension(input(%s, %%upload%%, 0,"
                          "'%s'), %s), %s)\n" %
                          (raw_array_schema, binary_tuple_schema,
                           sample_array_name, sample_array_name))
    else:
        script_file.write('#!/bin/bash\n\n')
        script_file.write('# Exit on error\nset -e\n\n')
        script_file.write('# Set the proper PATH\n')
        script_file.write("PATH=%s:$PATH\n\n" % opts.bindir)
        script_file.write("iquery -p %d -q \"create array %s%s\"\n" %
                          (opts.port, sample_array_name, sample_array_schema))
        script_file.write("iquery -p %d -an -q \"insert(redimension(input(%s, '%s'"
                          ", 0, '%s'), %s), %s)\"\n" %
                          (opts.port, raw_array_schema, out_file_name, binary_tuple_schema,
                           sample_array_name, sample_array_name))
