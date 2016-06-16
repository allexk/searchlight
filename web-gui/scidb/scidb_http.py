#!/usr/bin/python
"""This module contains function for accessing SciDB via its HTTP API.

All actions go throught the SciDBConnection class, which encapsulates a single SciDB session. See the
description of class for more detail.

Limitations: Currently, no HTTPS is supported. The main reason is that the current deployment does not use it.
It can be easily added later by request.
"""

__title_ = 'scidb_http'
__author__ = 'Alexander Kalinin <akalinin@cs.brown.edu>'
__version__ = '1.0.0'

import os
import sys
import subprocess
import requests
from requests.auth import HTTPDigestAuth
from requests.packages.urllib3.util import Retry
from scidb4py import Connection as SciDBBinaryConnection
from scidb4py import InternalError as SciDBBinaryInternalError


class SciDBConnection(object):
    """Class encapsulating a single SciDB session.

    A SciDB session consists of a series of requests (e.g., queries, file uploads, etc.) defined by the
    methods. All communication with SciDB is made via its HTTP API.

    All communications is handled by Python with one notable exception. File uploads are handled by an
    external call to curl. requests might support efficient file streaming, in which case uploading will be
    transferred to requests later.
    """
    def __init__(self, host, port, user, password):
        """Create a new SciDB session for executing queries.

        No connection is made at this point. User must call connect() for that.

        Params:
            host -- host to connect to (full DNS name)
            port -- port for connection
            user -- user name
            password -- password in clear text
        """
        # connection params
        self._user = user
        self._password = password
        self._host = host
        self._base_url = 'http://%s:%d' % (host, port)
        self._session = requests.Session()
        self._session.auth = HTTPDigestAuth(user, password)
        self._binary_connection = None

    def session_id(self):
        """Return SciDB session id.

        If no session is open, return None.
        """
        if 'id' in self._session.params:
            return self._session.params['id']
        else:
            return None

    def connect(self):
        """Connect to SciDB and create a new session."""
        print 'Connecting to %s...' % self._base_url
        body = self._get('/new_session', {})
        self._session.params['id'] = body.strip()
        print 'Success, session_id=%s' % self.session_id()

    def query(self, query_str, need_result=True, res_format='dcsv'):
        """Execute a SciDB query.

        A session must be already opened by connect(). Note, this function does not return th result.
        It must be retrieved later by get_result().

        Params:
            query_str -- SciDB query string
            need_result -- True, if the result is needed later; False, otherwise.
            res_format -- format of the result
        Return:
            query_id -- string
        """
        if not self.session_id():
            raise RuntimeError('Cannot query SciDB without connecting first')
        if not need_result:
            save_format = None
        else:
            save_format = res_format
        print "Executing query '%s'..." % query_str
        body = self._get('/execute_query', {'query': query_str, 'save': save_format})
        query_id = body.strip()
        print 'Success, query_id=%s' % query_id
        return query_id

    def query_interactive(self, query_str):
        """Execute SciDB query with interactive result fetching.

        The main difference between this function and query() is the mode
        of connection and result retrieval. This function creates an additional
        connection directly to SciDB (i.e., as iquery would do) and retrieves
        results by chunks.

        Returns:
            result array that can be iterated over(see scidb4py)
        """
        self._binary_connection = SciDBBinaryConnection(self._host)
        print "Connecting to SciDB via binary on port 1239..."
        self._binary_connection.open()
        print "Executing query '%s' via the binary connection..." % query_str
        res_array = self._binary_connection.execute(query_str, afl=True)
        return res_array

    def cancel_query(self):
        """Cancels current query associated either with this connection.

        This function will first try to cancel the binary connection query. Then, it will cancel
        the HTTP connection (session) query as well.
        """
        if self._binary_connection:
            try:
                print 'Cancelling binary connection query...'
                self._binary_connection.cancel()
            except SciDBBinaryInternalError:
                # No active query probably
                pass
        if self.session_id():
            print 'Cancelling session query...'
            self._get('/cancel', {})
        print 'Done.'

    def get_result(self, query_id):
        """Retrieve result for a previously executed query.

        For now the result is restricted to the last query executed. It seems SciDB HTTP API does not support
        result retrieval for an arbitrary query (e.g., by its id).

        Params:
            query_id -- query id, from query()
        Return:
            the entire result in a single string
        """
        if not self.session_id():
            raise RuntimeError('Cannot query SciDB without connecting first')
        # Actually, it seems query_id is not required for the result, but it enforces the user to run query() first
        print 'Retrieving result for query %s...' % query_id
        body = self._get('/read_lines', {'n': '0'})  # n means the number of lines (0 -- all lines)
        print 'Success'
        return body

    def close(self):
        """Close current session.

        If no session is open, do nothing.
        """
        if self._binary_connection:
            print 'Closing SciDB binary connection...'
            self._binary_connection.close()
            self._binary_connection = None
        if self.session_id():
            print 'Closing SciDB session %s...' % self.session_id()
            try:
                self._get('/release_session', {})
            except requests.HTTPError:
                # for some errors SciDB cleans up the session itself
                pass
            self._session.params.pop('id')
        print 'Success'

    def upload(self, file_path):
        """Upload a file to SciDB.

        This class uses an external curl process to transfer the file due to efficiency reasons. Note,
        SciDB returns a new, internal, file name. This file name must be used in subsequent queries.

        Params:
            file_path -- path to the file to upload
        Return:
            SciDB internal file name of the uploaded file
        """
        if not self.session_id():
            return
        if not os.path.isfile(file_path):
            raise ValueError('%s must specify a valid file' % file_path)

        _, file_name = os.path.split(file_path)
        print "Uploading '%s' to SciDB..." % file_path
        url = '%s/upload_file?id=%s' % (self._base_url, self.session_id())
        curl_args = ['curl', '--digest', '-u', '%s:%s' % (self._user, self._password),
                     '-F', '"file=@%s;filename=%s"' % (file_path, file_name), url]
        scidb_file_name = subprocess.check_output(curl_args).strip()
        print 'Success, scidb_file_name=%s' % scidb_file_name
        return scidb_file_name

    def _get(self, resource, url_params, retries=5):
        """Make an HTTP request and handle errors, if any."""
        url = self._base_url + resource
        try:
            r = self._session.get(url, params=url_params)
        except requests.exceptions.ConnectionError:
            # try to reconnect in case the connection is dropped
            if retries > 0:
                return self._get(resource, url_params, retries - 1)
        body_text = r.text
        if r.status_code != 200:
            error_string = '%d %s' % (r.status_code, r.reason)
            if body_text:
                # SciDB passes syntax errors in the body
                error_string += '\n%s' % body_text
            raise requests.HTTPError(error_string)
        return body_text

    def __del__(self):
        """Close the connection on deletion."""
        self.close()


def _main():
    """Runs if the module is executed directly as a script."""
    # command line parser
    import argparse
    parser = argparse.ArgumentParser(
        description="Connects to SciDB and runs commands according to the script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # command line parameters
    parser.add_argument('--port', metavar='N', type=int, default=8080, help='Port to connect to')
    parser.add_argument('--user', metavar='name', default='SciDBUser', help='User to connect as')
    parser.add_argument('--password', metavar='path', help='Path to a file with the password')
    parser.add_argument('--db-path', metavar='path', default='/home/gridsan/groups/databases',
                        help='Path to the databases tree')
    parser.add_argument('host', help='Host to connect to or the name of the instance')
    parser.add_argument('script', help='Script containing the commands to execute')
    # parse
    opts = parser.parse_args()

    # figure our the host
    if '.' in opts.host or opts.host == 'localhost':
        host = opts.host
    else:
        try:
            with open(os.path.join(opts.db_path, opts.host, 'dnsname'), 'r') as dns_file:
                host = dns_file.read().strip()
        except IOError:
            print 'Cannot figure out the host for %s. ' % opts.host
            return 1
    # figure out the password
    password_path = opts.password
    if not password_path:
        scidb_name = opts.host
        if '.' in scidb_name:
            scidb_name = scidb_name.split('.', 1)[0]
        password_path = os.path.join(opts.db_path, scidb_name, 'scidb_shim_password.txt')
    # read the password
    with open(password_path, 'rb') as password_file:
        password = password_file.read().strip()

    # connect
    scidb = SciDBConnection(host, opts.port, opts.user, password)
    try:
        scidb.connect()
        with open(opts.script, 'r') as script_file:
            query_ids = []
            uploaded_files = []
            for command in script_file:
                command = command.strip()
                if command.startswith('#'):  # script comment
                    continue
                action, command = command.split(':', 1)
                if action == 'query' or action == 'query_res':
                    if action == 'query':
                        need_result = False
                    else:
                        need_result = True
                    if '%upload%' in command:
                        command = command.replace('%upload%', "%s" % uploaded_files[-1])
                    query_id = scidb.query(command, need_result)
                    query_ids.append(query_id)
                elif action == 'query_int':
                    if '%upload%' in command:
                        command = command.replace('%upload%', "%s" % uploaded_files[-1])
                    res_array = scidb.query_interactive(command)
                    for pos, val in res_array:
                        # in this example assume single coordinate/attribute
                        print '%d - %s' % (pos[0], str(val[0]))
                elif action == 'upload':
                    file_name = scidb.upload(command)
                    uploaded_files.append(file_name)
                elif action == 'result':
                    if not query_ids:
                        raise ValueError('No query to get results for!')
                    else:
                        query_id = query_ids[-1]
                    with open(command, 'wb') as result_file:
                        result_file.write(scidb.get_result(query_id))
                        print 'Saved result for query %s into %s' % (query_id, command)
                else:
                    raise ValueError('Unknown action: %s' % action)
    finally:
        scidb.close()
    return 0

if __name__ == '__main__':
    sys.exit(_main())
