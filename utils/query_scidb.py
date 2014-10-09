#!/usr/bin/python
from __future__ import print_function

import os
import sys
import traceback
import datetime

# impoer scidb
sys.path.append('/opt/scidb/14.3/lib')
import scidbapi as scidb

def handle_exception(inst, exit, op=None):
    "Prints exception info and terminates if needed"
    traceback.print_exc()
    if op:
        print("Exception while", op, file=sys.stderr)
    print("     Exception Type:", type(inst), file=sys.stderr)
    print("     Exception Value: %r" % inst, end='\n\n', file=sys.stderr)
    if exit:
        sys.exit(2)

def get_query(query_file):
    "Reads the query from the specified file"
    with open(query_file, 'r') as q:
        query = q.read().strip()
    return query

def print_stats(times):
    "Print times statistics about the query execution"
    print("\nQuery started at:", times[0])
    print("Query finished at:", times[-1])
    print("Query total time:", times[-1] - times[0])
    
    if len(times) > 2:
        chunk_durs = [times[i] - times[i - 1] for i in range(1, len(times) - 1)]
            
        print("Got %d result chunks..." % len(chunk_durs))
        print("Minimum result delay:", min(chunk_durs))
        print("Maximum result delay:", max(chunk_durs))
        print("Average result delay:",
              sum(chunk_durs, datetime.timedelta()) / len(chunk_durs))
    else:
        print("No results found --- no duration stats!")

def execute_query(query):
    "Executes the query and times start, end and result times"
    times = []
    try:
        db = scidb.connect("localhost", 1239)
    except Exception, inst:
        # terminate
        handle_exception(inst, True, op="connecting")

    print("Executing:", query)
    try:
        times.append(datetime.datetime.now())
        result = db.executeQuery(query, "afl")
    except Exception, inst:
        # terminate
        handle_exception(inst, True, op="Performing select")

    # More precisely, we measure chunk times, since all result in a chunk
    # can be output with minimal delays. This should be enough to measure
    # interactive performance.
    desc = result.array.getArrayDesc()
    dims = desc.getDimensions()
    attrs = desc.getAttributes()

    # print result's schema
    print("Dimensions:")
    for i in range(dims.size()):
        print("\tDimension[%d] = %d:%d,%d,%d\n" % \
            (i, dims[i].getStartMin(), dims[i].getEndMax(),
             dims[i].getChunkInterval(), dims[i].getChunkOverlap()))
    print("Attributes:")
    for i in range(attrs.size()):
        print("\tAttribute %d %s %s" % \
            (attrs[i].getId(), attrs[i].getName(), attrs[i].getType()))

    # attribute partitioning: an iterator per attribute
    iters = []
    for i in range(attrs.size()):
        if not attrs[i].getName() == "EmptyTag":
            attrid = attrs[i].getId()
            iters.append(result.array.getConstIterator(attrid))

    # retrieve chunks and time them
    while not iters[0].end():
        # get and time another portion of the result
        chunks = []
        for i in range(len(iters)):
            chunks.append(iters[i].getChunk())
        times.append(datetime.datetime.now())

        for i in range(len(iters)):
            iters[i].increment_to_next()

    db.completeQuery(result.queryID)
    times.append(datetime.datetime.now())
    db.disconnect()  # Disconnect from the SciDB server.
    return times

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python query_scidb.py <query_file>")
        sys.exit(1)

    query = get_query(sys.argv[1])
    times = execute_query(query)
    print_stats(times)
    sys.exit(0)
