__author__ = 'akalinin'

import json
import os
import sys
import stat
import tempfile
from scidb.scidb_http import SciDBConnection
from ConfigParser import SafeConfigParser, NoOptionError
import shutil
import threading
import subprocess

# Hackish way to use pickled mimic index ("unused" imports are needed for unpickling!)
from mimic import MIMICWaveformData, Patient, PatientRecord
import mimic_index
MIMIC_CACHE = mimic_index.MIMICIndex("mimic.cache")

# Flask application
import flask
app = flask.Flask("searchlight")
app.config.from_pyfile("flask.cfg")

# create a logger
import logging
logger = logging.getLogger('searchlight')
logger_stream_handler = logging.StreamHandler()
logger_stream_handler.setFormatter(
    logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
logger.addHandler(logger_stream_handler)
logger.setLevel(logging.INFO)
# for logging all events (e.g., from modules) uncomment below
# logging.getLogger().addHandler(logger_stream_handler)
# logging.getLogger().setLevel(logging.DEBUG)
del logger_stream_handler

# to debug requests/responses
# try:
#     import http.client as http_client
# except ImportError:
#     # Python 2
#     import httplib as http_client
# http_client.HTTPConnection.debuglevel = 1

# Some default params (TODO: put them as params)
TASK_ARRAY = "mimic"
SIM_SYNOPSIS_PREFIX = "paa_128"

AVG_SYNOPSIS_RESOLUTIONS = ['1x1000', '1x100']
SIM_SYNOPSIS_RESOLUTIONS = ['1x200']

AVG_TASK_LIB = 'searchlight_mimic'
SIM_TASK_LIB = 'searchlight_dist'

AVG_TASK_FUN = 'MimicAvg'
SIM_TASK_FUN = 'Dist'


class SearchlightError(RuntimeError):
    """Searchlight exception."""
    pass


class SciDBConnectionConfig(object):
    def __init__(self, config_path):
        self._read_config(config_path)
        self.scidb = SciDBConnection(self.host, self.http_port, self.user,
                                     self.password)
        super(SciDBConnectionConfig, self).__init__()

    def _read_config(self, config_path):
        config = SafeConfigParser()
        config.read(config_path)

        # standard params
        self.host = config.get('scidb', 'host')
        self.user = config.get('scidb', 'user')
        self.http_port = int(config.get('scidb', 'http_port'))
        self.scidb_port = int(config.get('scidb', 'scidb_port'))
        self.limit = int(config.get('scidb', 'limit'))
        try:
            self.shared_dir = config.get('scidb', 'shared_dir')
            # assume list of nodes
            self.shared_nodes = eval(config.get('scidb', 'shared_nodes'))
        except NoOptionError:
            self.shared_dir = None
            self.shared_nodes = []

        # password we get from the specified file
        password_path = config.get('scidb', 'passwd_file')
        try:
            with open(password_path, 'rb') as password_file:
                self.password = password_file.read().strip()
        except IOError:
            logger.error('SciDB password file not found: ' + password_path)


class SciDBQueryOnline(SciDBConnectionConfig):
    def __init__(self, config_path):
        super(SciDBQueryOnline, self).__init__(config_path)
        self.logger = logging.getLogger('searchlight.scidb_query')
        self.scidb_query_params = None
        self.query_array = None
        self.cancelled = False
        self.error = None
        self.results = []

    def prepare_query(self, array, sample_arrays, task_name, task_fun,
                      task_params_file):
        self.scidb_query_params = {'array': array,
                                   'sample_arrays': sample_arrays,
                                   'task_name': task_name,
                                   'task_fun': task_fun,
                                   'task_params_file': task_params_file
                                   }

    def upload_file(self, filename):
        self.logger.info('Uploading file %s to all nodes...' % filename)
        if self.shared_dir:
            # copy locally
            shutil.copy(filename, self.shared_dir)
            # and remotely
            if len(self.shared_nodes):
                scp_args = ['scp', filename]
                for node in self.shared_nodes:
                    self.logger.info("copying file to '%s'..." % node)
                    print subprocess.check_output(scp_args + ['%s:%s' % (node, self.shared_dir)],
                                                  stderr=subprocess.STDOUT)
            new_filename = os.path.join(self.shared_dir, os.path.basename(filename))
        else:
            self.logger.info('connecting to %s:%d' % (self.host, self.http_port))
            # We need to connect only if we need to upload files via HTTP
            self.scidb.connect()
            # remote copy (curl via scidb shim)
            new_filename = self.scidb.upload(filename)
        self.logger.info('upload finished: ' + new_filename)
        return new_filename

    def start_query(self):
        if not self.scidb_query_params:
            self.logger.error('run_query() is called before prepare_query()')
            self.error = RuntimeError("query not inited")
            return
        try:
            # upload the query file
            task_params_file = self.scidb_query_params['task_params_file']
            scidb_task_file = self.upload_file(task_params_file)
            # create the query
            sample_arrays_str = ', '.join(
                ['depart(%s)' % x
                 for x in self.scidb_query_params['sample_arrays']])
            query_str = "searchlight(depart(%s), %s, '%s:%s:%s')" %\
                        (self.scidb_query_params['array'],
                         sample_arrays_str,
                         self.scidb_query_params['task_name'],
                         self.scidb_query_params['task_fun'],
                         scidb_task_file)
            # query (via scidb4py)
            self.logger.info('querying with: ' + query_str)
            self.query_array = iter(self.scidb.query_interactive(query_str))
        except:
            self.logger.error('could not start query: ' + str(sys.exc_info()))
            self.error = sys.exc_info()

    def next_result(self):
        if not self.query_array:
            self.logger.info('no current array, possible run_query() error...')
            self._finish_query()
            self.error = RuntimeError("query not started")
            return None
        try:
            if self.cancelled or len(self.results) >= self.limit:
                self.logger.info('query cancelled by user or limit reached (%d)' % self.limit)
                try:
                    self.scidb.cancel_query()
                except:
                    pass
                self._finish_query()
            else:
                (pos, val) = self.query_array.next()
                # assume the only attribute -- the result string
                res = self._parse_str_result(str(val[0]))
                self.results.append(res)
                return res
        except StopIteration:
            # close the connection
            self.logger.info('no more results: query finished')
            self._finish_query()
        return None

    def cancel_query(self):
        self.cancelled = True

    def _finish_query(self):
        if self.query_array:
            self.logger.info('closing SciDB connection')
            try:
                self.scidb.close()
            except:
                pass

    def get_result(self, ind):
        if ind < 0 or ind > len(self.results):
            raise ValueError("result index out of bounds")
        return self.results[ind]

    @staticmethod
    def _parse_str_result(res_str):
        """Parses result into a dict.

        A waveform result is traditionally (id, start, len)
        """
        res = {}
        for r in res_str.strip().split(', '):
            (name, val) = r.split('=')
            res[name] = int(val)
        return res

    def __del__(self):
        self._finish_query()


class MockSciDBQueryOnline(SciDBQueryOnline):
    def __init__(self, config_path):
        super(MockSciDBQueryOnline, self).__init__(config_path)
        self.logger = logging.getLogger('searchlight.mock_scidb_query')

    def start_query(self):
        if not self.scidb_query_params:
            self.logger.error('run_query() is called before prepare_query()')
            self.error = RuntimeError("query not inited")
            return
        try:
            # upload the query file
            task_params_file = self.scidb_query_params['task_params_file']
            self.logger.info('uploading task file ' + task_params_file)
            if self.shared_dir:
                # copy locally
                shutil.copy(task_params_file, self.shared_dir)
                scidb_task_file = os.path.join(self.shared_dir,
                                               os.path.basename(
                                                   task_params_file))
            else:
                # connect
                self.logger.info('connecting to %s:%d' % (self.host, self.http_port))
                # remote copy (curl via scidb shim)
                scidb_task_file = "<remote_file>"
            self.logger.info('upload finished: ' + scidb_task_file)
            # create the query
            sample_arrays_str = ', '.join(
                ['depart(%s)' % x
                 for x in self.scidb_query_params['sample_arrays']])
            query_str = "searchlight(depart(%s), %s, '%s:%s:%s')" %\
                        (self.scidb_query_params['array'],
                         sample_arrays_str,
                         self.scidb_query_params['task_name'],
                         self.scidb_query_params['task_fun'],
                         scidb_task_file)
            # query (via scidb4py)
            self.logger.info('querying with: ' + query_str)
            self.query_array = iter(["id=0, time=1, len=5",
                                     "id=1, time=10, len=50",
                                     "id=0, time=20, len=55"])
        except:
            self.logger.error('could not start query: ' + str(sys.exc_info()))
            self.error = sys.exc_info()

    def next_result(self):
        if not self.query_array:
            self.logger.info('no current array, possible run_query() error...')
            self._finish_query()
            self.error = RuntimeError("query not started")
            return None
        try:
            if self.cancelled:
                self.logger.info('query cancelled by user')
                self._finish_query()
            else:
                val = self.query_array.next()
                # assume the only attribute -- the result string
                res = self._parse_str_result(str(val))
                self.results.append(res)
                return res
        except StopIteration:
            # close the connection
            self.logger.info('no more results: query finished')
            self._finish_query()
        return None

    def cancel_query(self):
        self.cancelled = True

    def _finish_query(self):
        if self.query_array:
            self.logger.info('closing SciDB connection')

    def get_result(self, ind):
        if ind < 0 or ind > len(self.results):
            raise ValueError("result index out of bounds")
        return self.results[ind]

    @staticmethod
    def _parse_str_result(res_str):
        """Parses result into a dict.

        A waveform result is traditionally (id, start, len)
        """
        res = {}
        for r in res_str.strip().split(', '):
            (name, val) = r.split('=')
            res[name] = int(val)
        return res

    def __del__(self):
        self._finish_query()


class QueryJSONParamsHandler(object):
    """
    Class to handle getting/serializing query params via JSON.

    When working with params this class assumes path-like syntax, e.g.,
    mimic.sl.db. All paths are directed to the mimic JSON element.
    For example, a.b will retrieve mimic.a.b.
    """
    def __init__(self, f):
        self.query_json = json.load(f)

    def set_params(self, params):
        for (param, val) in params.iteritems():
            param_path = param.split('.')
            param_ref = self.get_param(param_path[:-1])
            param_ref[param_path[-1]] = val

    def get_param(self, path):
        if isinstance(path, str):
            path = path.split('.')
        #  param_val = self.query_json['mimic']
        param_val = self.query_json
        for step in path:
            param_val = param_val[step]
        return param_val

    def dump_to_temp(self):
        (fd, name) = tempfile.mkstemp(dir="tmp")
        os.fchmod(fd, stat.S_IRUSR | stat.S_IWUSR | stat.S_IROTH)
        with os.fdopen(fd, 'wb') as f:
            json.dump(self.query_json, f, sort_keys=True, indent=4,
                      separators=(',', ': '))
        return name

    def __str__(self):
        """String representation."""
        return str(self.query_json)


class SciDBWaveformQuery(SciDBConnectionConfig):
    def __init__(self, config_path):
        super(SciDBWaveformQuery, self).__init__(config_path)
        self.logger = logging.getLogger('searchlight.scidb_waveform')

    def get_waveform(self, array, signal, record_id, start, length):
        self.logger.info('Connecting to SciDB (waveform)...')
        res = []
        try:
            query = 'between(%s, %d, %d, %d, %d)' % (array, record_id, start,
                                                     record_id,
                                                     start + length - 1)
            self.logger.info('Retrieving waveform: ' + query)
            res_array = self.scidb.query_interactive(query)
            current_pos = start
            for (coords, attrs) in res_array:
                pos = int(coords['tick'])  # the timeline coordinate
                for i in range(current_pos, pos):
                    # empty elements -- assume 0
                    res.append(0.0)
                current_pos = pos + 1
                res.append(float(attrs[signal]))
            self.logger.info('Closing waveform SciDB connection...')
            self.scidb.close()
        except:
            # Cannot do anything about it. Just return empty waveform
            self.logger.error(str(sys.exc_info()))
        return res


class MockSciDBWaveformQuery(SciDBConnectionConfig):
    def get_waveform(self, array, signal, record_id, start, length):
        import random
        return [random.randint(50, 100) for _ in range(length)]


class QueryHandler(object):
    def __init__(self):
        self._lock = threading.Lock()
        self._id = 0
        self._queries = {}

    def set_query(self, query):
        with self._lock:
            qid = str(self._id)
            self._id += 1
            self._queries[qid] = query
        logger.info("query registered: " + str(qid))
        return qid

    def get_query(self, query_id):
        with self._lock:
            if query_id not in self._queries:
                raise SearchlightError("no query: " + str(query_id))
            return self._queries[query_id]

    def remove_query(self, qid):
        with self._lock:
            if qid in self._queries:
                logger.info("removed query: " + str(qid))
                del self._queries[qid]
            else:
                logger.warning("no query with id: " + str(qid))
global_query_handler = QueryHandler()


def cleanup_query(query_id):
    """Remove query from the list of active queries."""
    logger.info("cleaning up query: " + query_id)
    global_query_handler.remove_query(query_id)


def parse_json_request(request):
    """Parse JSON request.

    Params:
        :param request server request
    Return:
        tuple: (query_id, query/None)
    """
    # parse JSON request
    data_json = request.get_json()
    if data_json is None:
        raise SearchlightError("Expecting JSON query")
    # The query has two parts: "query"...
    try:
        query_json = data_json["query"]
        if query_json:
            json_try_int(query_json)
        logger.info("Got query: " + str(query_json))
    except KeyError:
        query_json = None
    # and "id"
    try:
        query_id = data_json["query_id"]
        logger.info("Got query_id: " + str(query_id))
    except KeyError:
        # No query_id field or null
        query_id = None
    # and aux_data
    try:
        query_data = data_json["aux_data"]
        if query_data:
            logger.info("Got aux_data of length %d" % len(query_data))
    except KeyError:
        # No aux data for the query or null
        query_data = None
    return query_id, query_json, query_data


def json_try_int(js):
    """Try to go through JSON and convert values to ints."""
    for key in js.keys():
        try:
            val_i = int(js[key])
            js[key] = val_i
        except ValueError:
            pass


def json_transform(js, keep_keys, prefix_rename):
    """Go through the JSON, keep some keys and rename prefixes"""
    new_js = {}
    for key in js.keys():
        if keep_keys is None or key in keep_keys:
            new_key = key
            if prefix_rename is not None and key.startswith(prefix_rename[0]):
                new_key = key.replace(prefix_rename[0], prefix_rename[1], 1)
            new_js[new_key] = js[key]
    return new_js


@app.errorhandler(SearchlightError)
def respond_with_error(err):
    response = flask.jsonify({"status": "error", "message": str(err)})
    response.status_code = 400
    logger.error("error occured:" + str(err))
    return response


@app.route('/')
def root():
    return app.send_static_file('searchlight.html')


@app.route('/js/<path:filepath>')
def send_js(filepath):
    return flask.send_from_directory('static/js', filepath)


@app.route('/css/<path:filepath>')
def send_css(filepath):
    return flask.send_from_directory('static/css', filepath)


@app.route('/fonts/<path:filepath>')
def send_fonts(filepath):
    return flask.send_from_directory('static/fonts', filepath)


@app.route("/query", methods=['POST'])
def start_query():
    # parse JSON request
    query_json = parse_json_request(flask.request)[1]

    with open("mimic_tmpl.json", "rb") as f:
        query_handler = QueryJSONParamsHandler(f)
    query_handler.set_params(query_json)
    query_file = query_handler.dump_to_temp()

    # create SciDB query
    sl_query = SciDBQueryOnline("config.ini")
    synopsis_arrays = ['_'.join([TASK_ARRAY,
                                 query_handler.get_param("mimic.signal"),
                                 syn_res]) for syn_res in AVG_SYNOPSIS_RESOLUTIONS]
    sl_query.prepare_query(TASK_ARRAY, synopsis_arrays, AVG_TASK_LIB, AVG_TASK_FUN,
                           query_file)
    sl_query.start_query()

    # response
    if not sl_query.error:
        query_id = global_query_handler.set_query(sl_query)
        return flask.jsonify({"status": "ok", "query_id": query_id})
    else:
        raise SearchlightError("Exception when running the query")


@app.route("/sim", methods=['POST'])
def start_sim_query():
    # parse JSON request
    _, query_json, query_data = parse_json_request(flask.request)
    if query_json is None or query_data is None:
        raise SearchlightError("Invalid sim request")
    # transform for the dist query
    dist_keys = {"mimic.l_id", "mimic.u_id", "mimic.l_time", "mimic.u_time", "mimic.dist", "mimic.step_time",
                 "mimic.signal", "relax.on", "relax.spec"}
    query_json = json_transform(query_json, dist_keys, ("mimic.", "dist."))
    logger.debug("Transformed query: %s" % str(query_json))
    # create SciDB query
    sl_query = SciDBQueryOnline("config.ini")
    # upload the sequence
    (data_fd, data_name) = tempfile.mkstemp(dir="tmp")
    os.fchmod(data_fd, stat.S_IRUSR | stat.S_IWUSR | stat.S_IROTH)
    with os.fdopen(data_fd, 'wb') as f:
        f.write(str(len(query_data)))
        for x in query_data:
            f.write(" %f" % x)
    query_json["dist.query"] = sl_query.upload_file(data_name)
    # prepare query
    with open("dist_tmpl.json", "rb") as f:
        query_handler = QueryJSONParamsHandler(f)
    query_handler.set_params(query_json)
    query_file = query_handler.dump_to_temp()
    logger.debug("Final query: %s" % str(query_handler))
    # create SciDB query
    synopsis_arrays = ['_'.join([SIM_SYNOPSIS_PREFIX, TASK_ARRAY,
                                 query_handler.get_param("dist.signal"),
                                 syn_res]) for syn_res in SIM_SYNOPSIS_RESOLUTIONS]
    logger.debug("Synopsis arrays will use: %s" % str(synopsis_arrays))
    sl_query.prepare_query(TASK_ARRAY, synopsis_arrays, SIM_TASK_LIB, SIM_TASK_FUN,
                           query_file)
    sl_query.start_query()
    # response
    if not sl_query.error:
        query_id = global_query_handler.set_query(sl_query)
        return flask.jsonify({"status": "ok", "query_id": query_id})
    else:
        raise SearchlightError("Exception when running the query")


# POST method since BigDawg requires it
@app.route("/result", methods=['POST'])
def next_result():
    query_id = parse_json_request(flask.request)[0]
    if query_id is None:
        raise SearchlightError("query_id is not specified")

    # get the next result (might block for a while)
    sl_query = global_query_handler.get_query(query_id)
    res_dict = sl_query.next_result()
    if not res_dict:
        cleanup_query(query_id)
        if sl_query.error:
            raise SearchlightError("Unexpected query error")
        else:
            res_dict = {"eof": True}
    else:
        # now we should determine the real time and subject id
        record_id = int(res_dict["id"])
        start_tick = int(res_dict["time"])
        patient_id, record_name, pretty_time = MIMIC_CACHE.find_segment(record_id, start_tick)
        patient_id = int(patient_id[1:])  # assuming sxxx format (e.g., s00124)
        res_dict["sid"] = patient_id
        res_dict["pretty_time"] = pretty_time
        res_dict = {"status": "ok", "result": res_dict, "eof": False}
    return flask.jsonify(res_dict)


# POST method since BigDawg requires it
@app.route("/waveform", methods=["POST"])
def waveform():
    # parse request
    query_json = parse_json_request(flask.request)[1]
    for p in ["signal", "id", "time", "len"]:
        if p not in query_json:
            raise SearchlightError("missing parameter " + p)
    logger.debug("Got request: " + str(query_json))

    # Get waveform
    waveform_getter = SciDBWaveformQuery("config.ini")
    res = waveform_getter.get_waveform(TASK_ARRAY, query_json["signal"],
                                       query_json["id"], query_json["time"],
                                       query_json["len"])
    if len(res) > 0:
        return flask.jsonify({"status": "ok", "waveform": res})
    else:
        raise SearchlightError("cannot retrieve waveform")


@app.route("/cancel", methods=['POST'])
def cancel_query():
    query_id = parse_json_request(flask.request)[0]
    if query_id is None:
        raise SearchlightError("query_id is not specified")

    # cancel and return immediately
    try:
        sl_query = global_query_handler.get_query(query_id)
        sl_query.cancel_query()
    except SearchlightError:
        # ignore, nothing fatal here
        pass
    return flask.jsonify({"status": "ok"})

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5000, threaded=True)
