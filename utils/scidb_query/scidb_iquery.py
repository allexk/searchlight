#!/usr/bin/python
"""This module contains function for accessing SciDB via its HTTP API.

All actions go throught the SciDBConnection class, which encapsulates a single SciDB session. See the
description of class for more detail.

Limitations: Currently, no HTTPS is supported. The main reason is that the current deployment does not use it.
It can be easily added later by request.
"""

__title_ = 'scidb_iquery'
__author__ = 'Alexander Kalinin <akalinin@cs.brown.edu>'
__version__ = '1.0.0'

import os
import subprocess
import cStringIO as StringIO


class SciDBConnection(object):
    """Class encapsulating a single SciDB session.

    A SciDB session consists of a series of requests (e.g., queries, file uploads, etc.) defined by the
    methods. All communication with SciDB is made via iquery utility.
    """
    def __init__(self, scidb_bin_path, scidb_port):
        """Create a new SciDB session for executing queries.

        Params:
            scidb_bin_path -- path to SciDB binaries
            scidb_port -- SciDB port to connect to
        """
        # connection params
        self._iquery = os.path.join(scidb_bin_path, 'iquery')
        if not os.path.isfile(self._iquery):
            raise ValueError("Cannot find iquery binary in %s" % scidb_bin_path)
        self._query_results = {}
        self._iquery_port = scidb_port

    def session_id(self):
        """Return SciDB session id.

        If no session is open, return None.
        """
        # For iquery it's always the same "session"
        return 0

    def connect(self):
        """Connect to SciDB and create a new session.

        Return:
            session_id -- string
        """
        pass

    def query(self, query_str, need_result=True, res_format='dcsv'):
        """Execute a SciDB query.

        A session must be already opened by connect(). Note, this function does not return th result.
        It must be retrieved later by get_result().

        Params:
            query_str -- SciDB query string
            need_result -- True, if the result is needed later; False, otherwise
            res_format -- format of the result
        Return:
            query_id -- string
        """
        iquery_args = [self._iquery, '-p', str(self._iquery_port), '-a']
        if not need_result:
            iquery_args.append('-n')
        else:
            iquery_args.extend(['-o', res_format])
        iquery_args.extend(['-q', query_str])
        query_id = len(self._query_results)
        print "Executing query '%s'..." % query_str
        result = subprocess.check_output(iquery_args, stderr=subprocess.STDOUT)
        self._query_results[query_id] = result
        print 'Success, query_id=%s' % query_id
        return query_id

    def get_result(self, query_id):
        """Retrieve result for a previously executed query.

        Params:
            query_id -- query id, from query()
        Return:
            the entire result in a single string
        """
        if query_id not in self._query_results:
            raise ValueError('There is no result available for: %d' % query_id)
        print 'Retrieving result for query %s...' % query_id
        return self._query_results[query_id]

    def close(self):
        """Close current session.

        If no session is open, do nothing.
        """
        pass

    def upload(self, file_path):
        """Upload a file to SciDB.

        Params:
            file_path -- path to the file to upload
        Return:
            SciDB internal file name of the uploaded file
        """
        if not os.path.isfile(file_path):
            raise ValueError('%s must specify a valid file' % file_path)
        # for iquery just leave the file at the same place
        return file_path

