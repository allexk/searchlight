__author__ = 'akalinin'

from PyQt5 import QtWidgets, QtCore
from main_window import Ui_MainWindow

import json
import os
import sys
import tempfile
from scidb_http import SciDBConnection
from ConfigParser import SafeConfigParser
import logging

# Some default params (TODO: put them as params)
SYNOPSIS_RESOLUTIONS = ['1x1000', '1x100']
#SYNOPSIS_RESOLUTIONS = ['1x1000']
TASK_LIB = 'searchlight_mimic'
TASK_FUN = 'MimicAvg'

# create a logger
logger = logging.getLogger('searchlight')
logger_stream_handler = logging.StreamHandler()
logger_stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
#logging.getLogger().addHandler(logger_stream_handler)
#logging.getLogger().setLevel(logging.DEBUG)
logger.addHandler(logger_stream_handler)
logger.setLevel(logging.INFO)
del logger_stream_handler

class SciDBConnectionConfig(object):
    def __init__(self, config_path):
        self._read_config(config_path)
        self.scidb = SciDBConnection(self.host, self.http_port, self.user, self.password)
        super(SciDBConnectionConfig, self).__init__()

    def _read_config(self, config_path):
        config = SafeConfigParser()
        config.read(config_path)

        # standard params
        self.host = config.get('scidb', 'host')
        self.user = config.get('scidb', 'user')
        self.http_port = int(config.get('scidb', 'http_port'))
        self.scidb_port = int(config.get('scidb', 'scidb_port'))

        # password we get from the specified file
        password_path = config.get('scidb', 'passwd_file')
        try:
            with open(password_path, 'rb') as password_file:
                self.password = password_file.read().strip()
        except IOError:
            logger.error('SciDB password file not found: ' + password_path)

class SciDBWaveformQuery(SciDBConnectionConfig):
    def __init__(self, config_path):
        super(SciDBWaveformQuery, self).__init__(config_path)
        self.logger = logging.getLogger('searchlight.scidb_waveform')

    def get_waveform(self, array, signal, record_id, start, length):
        self.logger.info('Connecting to SciDB (waveform)...')
        res = []
        try:
            query = 'between(%s, %d, %d, %d, %d)' % (array, record_id, start, record_id, start + length - 1)
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
            # probably cannot do anything about it. Worst -- return empty waveform
            self.logger.error(str(sys.exc_info()))
        return res

class WaveformResults(object):
    def __init__(self, config_path):
        # results is a dict: (id, start, len) --> waveform/None
        self.results = {}
        self.waveform_getter = SciDBWaveformQuery(config_path)

    def clear(self):
        self.results = {}

    def add_result(self, res_str):
        self.results[self._parse_str_result(res_str)] = None # no waveform

    def get_waveform(self, array, signal, res):
        if isinstance(res, str):
            res = self._parse_str_result(res)
        if res not in self.results:
            return None  # no such result
        if not self.results[res]:
            logger.info('retrieveing waveform for: ' + str(res))
            self.results[res] = self.waveform_getter.get_waveform(array, signal, *res)
        return self.results[res]

    @staticmethod
    def _parse_str_result(res_str):
        """Parses result into (id, start, len) tuple."""
        res = {}
        for r in res_str.strip().split(', '):
            (name, val) = r.split('=')
            res[name] = int(val)
        return res['id'], res['time'], res['len']

class SciDBQueryOnline(SciDBConnectionConfig, QtCore.QObject):
    def __init__(self, config_path):
        super(SciDBQueryOnline, self).__init__(config_path)
        self.logger = logging.getLogger('searchlight.scidb_query')
        self.scidb_query_params = None
        self.query_array = None
        self.cancelled = False

    def prepare_query(self, array, sample_arrays, task_name, task_fun, task_params_file):
        self.scidb_query_params = {'array': array,
                                   'sample_arrays': sample_arrays,
                                   'task_name': task_name,
                                   'task_fun': task_fun,
                                   'task_params_file': task_params_file
                                   }

    @QtCore.pyqtSlot()
    def run_query(self):
        if not self.scidb_query_params:
            self.logger.error('run_query() is called before prepare_query()')
            return
        try:
            # connect
            self.logger.info('connecting to %s:%d' % (self.host, self.http_port))
            self.scidb.connect()
            # upload the query file
            self.logger.info('uploading task file ' + self.scidb_query_params['task_params_file'])
            scidb_task_file = self.scidb.upload(self.scidb_query_params['task_params_file'])
            #scidb_task_file = '/home/scidb/sl_mimic.json'
            self.logger.info('upload finished: ' + scidb_task_file)
            # create the query
            sample_arrays_str = ', '.join(['depart(%s)' % x for x in self.scidb_query_params['sample_arrays']])
            query_str = "searchlight(depart(%s), %s, '%s:%s:%s')" % (self.scidb_query_params['array'],
                                                                     sample_arrays_str,
                                                                     self.scidb_query_params['task_name'],
                                                                     self.scidb_query_params['task_fun'],
                                                                     scidb_task_file)
            # query
            self.logger.info('querying with: ' + query_str)
            self.query_array = iter(self.scidb.query_interactive(query_str))
        except:
            self.logger.error('could not start query: ' + str(sys.exc_info()))

        self.scidb_query_params = None
        self.cancelled = False

    # custom signals
    result_ready = QtCore.pyqtSignal(str)
    query_finished = QtCore.pyqtSignal()

    @QtCore.pyqtSlot()
    def next_result(self):
        if not self.query_array:
            self.logger.info('no current array, possible run_query() error...')
            self._finish_query()
            return
        try:
            if self.cancelled:
                self.logger.info('query cancelled by user')
                try:
                    self.scidb.cancel_query()
                except:
                    pass
                self._finish_query()
            else:
                (pos, val) = self.query_array.next()
                self.result_ready.emit(str(val[0]))  # assume the only attribute -- the result string
        except StopIteration:
            # close the connection
            self.logger.info('no more results: query finished')
            self._finish_query()

    @QtCore.pyqtSlot()
    def cancel_query(self):
        self.cancelled = True

    def _finish_query(self):
        if self.query_array:
            self.logger.info('closing SciDB connection')
            try:
                self.scidb.close()
            except:
                pass
            self.query_array = None
        self.query_finished.emit()

class QueryJSONParamsHandler(object):
    """
    Class to handle getting/serializing query params via JSON.

    When working with params this class assumes path-like syntax, e.g., mimic.sl.db. All paths are
    directed to the mimic JSON element. For example, a.b will retrieve mimic.a.b.
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
        param_val = self.query_json['mimic']
        for step in path:
            param_val = param_val[step]
        return param_val

    def dump_to_temp(self):
        (fd, name) = tempfile.mkstemp(dir=os.getcwd())
        with os.fdopen(fd, 'wb') as f:
            json.dump(self.query_json, f, sort_keys=True, indent=4, separators=(',', ': '))
        return name

class Main(QtWidgets.QMainWindow, Ui_MainWindow):
    """Class responsible for the whole user interaction."""
    # A couple of signals to control the query thread
    _start_query = QtCore.pyqtSignal()
    _next_result = QtCore.pyqtSignal()
    _cancel_query = QtCore.pyqtSignal()

    class TmpFileHandler(object):
        """Stores path to the file with GC.

        When a new file is assigned, the previous one is unlinked. Also, it's
        unlinked on delete.
        """
        def __init__(self):
            self.file_path = None

        def __del__(self):
            if self.file_path:
                os.unlink(self.file_path)

        def add_file(self, file_path):
            if self.file_path:
                os.unlink(self.file_path)
            self.file_path = file_path

    # A couple of proxy classes to get/set values from/to widgets
    class ValueProxyInt(object):
        """Proxy to get the proper value from a value()-based widget."""
        def __init__(self, obj):
            self.obj = obj

        def value(self):
            """Get the value if the widget is enabled. 0 otherwise"""
            if self.obj.isEnabled():
                return self.obj.value()
            else:
                return 0

        def set_value(self, val):
            self.obj.setValue(val)

    class ComboBoxSetProxy(object):
        """Proxy to set the value to the ComboBox widget."""
        def __init__(self, obj):
            self.obj = obj

        def set(self, text):
            ind = self.obj.findText(text)
            if ind == -1:
                raise ValueError('Invalid value')
            self.obj.setCurrentIndex(ind)

    class TimeGetSetProxy(object):
        """Proxy to get/set the value from the LineEdit widget, treating it as time."""
        def __init__(self, getter, setter):
            self.getter = getter
            self.setter = setter

        def set(self, val):
            # for now just convert to str
            self.setter(str(val))

        def get(self):
            # for now just convert to str
            return str(self.getter())

    def __init__(self, config_path):
        super(Main, self).__init__()
        self.setupUi(self)

        self.query_handler = None  # to handle query params
        self.tmp_handler = self.TmpFileHandler()  # to handle temporary files
        self.params_dict = None  # query params dictionary
        self.in_query = False  # are we running the query?
        self.last_array = None  # last queried array
        self.last_signal = None  # last queried signal
        self.waveform_results = WaveformResults(config_path)  # to store results

        # Create the query thread
        self.scidb_query = SciDBQueryOnline(config_path)
        self.query_thread = QtCore.QThread()
        self.scidb_query.moveToThread(self.query_thread)
        self.query_thread.finished.connect(self.scidb_query.deleteLater)

        # Connect signals
        self.scidb_query.result_ready.connect(self._new_result)
        self.scidb_query.query_finished.connect(self._query_finished)
        self._cancel_query.connect(self.scidb_query.cancel_query)
        self._start_query.connect(self.scidb_query.run_query)
        self._next_result.connect(self.scidb_query.next_result)
        # Start (it'll run till we exit)
        self.query_thread.start()

        # Additional setup for params
        self.query_type.addItems(['mimic'])
        self.signal_type.addItems(['ABP'])

        # Various buttons
        self.nb_checkbox.stateChanged.connect(self.set_enable_neighborhood)
        self.start_button.clicked.connect(self.start_query)
        self.cancel_button.clicked.connect(self.cancel_query)

        # List widget waveform draw on double-click
        self.results_list.itemDoubleClicked.connect(self.draw_result)

        # Menu bar
        self.actionOpen.triggered.connect(self.open_file)
        self.actionExit.triggered.connect(QtWidgets.qApp.quit)

        # Some initialization
        self._set_dict()
        try:
            with open('mimic_tmpl.json', 'rb') as f:
                self.query_handler = QueryJSONParamsHandler(f)
                self._set_params()
        except (IOError, OSError):
            # Cannot open/parse that's fine -- no default params"""
            pass

    def __del__(self):
        self.query_thread.quit()
        self.query_thread.wait()

    def closeEvent(self, event):
        """Overrides closing event (the 'x' button) to avoid closing during query running."""
        if not self.in_query:
            super(Main, self).closeEvent(event)
        else:
            event.ignore()

    def _set_dict(self):
        """
        Creates the proper dictionary that contains information about params retrieval.

        The dict: original.query.param --> info tuple:
            (param type, getter, setter, error string).
        """
        self._val_proxies = {
            'neighborhood.l_size': self.ValueProxyInt(self.left_nb_size),
            'neighborhood.r_size': self.ValueProxyInt(self.right_nb_size),
            'neighborhood.left_max_diff': self.ValueProxyInt(self.left_nb_thr),
            'neighborhood.right_max_diff': self.ValueProxyInt(self.right_nb_thr),
            'signal': self.ComboBoxSetProxy(self.signal_type),
            'l_time': self.TimeGetSetProxy(self.left_search_bound.text, self.left_search_bound.setText),
            'r_time': self.TimeGetSetProxy(self.right_search_bound.text, self.right_search_bound.setText)
        }
        # param --> (type, getter, setter, pretty name)
        self.params_dict = {
            'l_time': (int, self._val_proxies['l_time'].get, self._val_proxies['l_time'].set, 'left search interval'),
            'u_time': (int, self._val_proxies['r_time'].get, self._val_proxies['r_time'].set, 'right search interval'),
            'avg_l': (int, self.min_avg.value, self.min_avg.setValue, 'min average'),
            'avg_u': (int, self.max_avg.value, self.max_avg.setValue, 'max average'),
            'len_l': (int, self.min_size.value, self.min_size.setValue, 'min length'),
            'len_u': (int, self.max_size.value, self.max_size.setValue, 'max length'),
            'step_time': (int, self.search_resolution.value, self.search_resolution.setValue, 'search resolution'),
            'signal': (str, self.signal_type.currentText, self._val_proxies['signal'].set, 'signal'),

            'neighborhood.l_size': (int, self._val_proxies['neighborhood.l_size'].value,
                                    self._val_proxies['neighborhood.l_size'].set_value, 'left neighborhood size'),
            'neighborhood.r_size': (int, self._val_proxies['neighborhood.r_size'].value,
                                    self._val_proxies['neighborhood.r_size'].set_value, 'right neighborhood size'),
            'neighborhood.left_max_diff': (int, self._val_proxies['neighborhood.left_max_diff'].value,
                                           self._val_proxies['neighborhood.left_max_diff'].set_value,
                                           'left maximum difference'),
            'neighborhood.right_max_diff': (int, self._val_proxies['neighborhood.right_max_diff'].value,
                                            self._val_proxies['neighborhood.right_max_diff'].set_value,
                                            'right maximum difference')
        }

    @QtCore.pyqtSlot()
    def open_file(self):
        """Retrieved a JSON file via an Open File Dialog"""
        (json_file, _) = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd(), 'JSON (*.json)')
        if json_file:
            try:
                with open(json_file, 'rb') as f:
                    self.query_handler = QueryJSONParamsHandler(f)
            except ValueError, e:
                QtWidgets.QMessageBox.warning(self, 'Error parsing JSON', str(e))
            self._set_params()

    @QtCore.pyqtSlot(int)
    def set_enable_neighborhood(self, state):
        """Enables/Disables setting of neighborhood params"""
        self.left_nb_size.setEnabled(state)
        self.right_nb_size.setEnabled(state)
        self.left_nb_thr.setEnabled(state)
        self.right_nb_thr.setEnabled(state)

    def _create_params_dict(self):
        """Creates a query parameter dict, getting the values from widgets."""
        res_dict = {}
        for (param, val_tuple) in self.params_dict.iteritems():
            (p_type, p_get, _, p_err_str) = val_tuple
            try:
                p_value = p_type(p_get())
            except ValueError, e:
                # Probably a wrong type in the box
                QtWidgets.QMessageBox.warning(self, 'Wrong parameter type',
                                              '%s: %s!' % (str(e), p_err_str))
                return None
            res_dict[param] = p_value
        return res_dict

    def _set_params(self):
        """Takes params from the query dict and sets widget values."""
        if not self.query_handler:
            return
        for (param, val_tuple) in self.params_dict.iteritems():
            (p_type, _, p_set, p_err_str) = val_tuple
            try:
                p_set(p_type(self.query_handler.get_param(param)))
            except ValueError, e:
                QtWidgets.QMessageBox.warning(self, 'Wrong parameter type',
                                              '%s: %s!' % (str(e), p_err_str))
            except KeyError, e:
                # No such param in JSON, pass
                pass

    @QtCore.pyqtSlot()
    def start_query(self):
        query_params = self._create_params_dict()
        if query_params:
            self.start_button.setDisabled(True)
            self.query_handler.set_params(query_params)
            tmp_name = self.query_handler.dump_to_temp()
            self.tmp_handler.add_file(tmp_name)
            self.last_array = str(self.query_type.currentText())
            self.last_signal = str(self.signal_type.currentText())
            synopsis_arrays = ['_'.join([self.last_array, self.last_signal, syn_res])
                               for syn_res in SYNOPSIS_RESOLUTIONS]
            self.scidb_query.prepare_query(self.last_array,
                                           synopsis_arrays,
                                           TASK_LIB,
                                           TASK_FUN,
                                           tmp_name)
            self.results_list.clear()
            self.waveform_results.clear()
            self._start_query.emit()
            self._next_result.emit()
            self.cancel_button.setEnabled(True)
            self.in_query = True

    @QtCore.pyqtSlot()
    def cancel_query(self):
        self._cancel_query.emit()
        self.cancel_button.setDisabled(True)

    @QtCore.pyqtSlot(str)
    def _new_result(self, res_str):
        """Displays the new result in the list."""
        self.results_list.addItem(res_str)
        self.waveform_results.add_result(res_str)
        self._next_result.emit()

    @QtCore.pyqtSlot()
    def _query_finished(self):
        """Called when the current query is finished."""
        self.start_button.setEnabled(True)
        self.cancel_button.setDisabled(True)
        self.in_query = False

    @QtCore.pyqtSlot(QtWidgets.QListWidgetItem)
    def draw_result(self, item):
        res_text = str(item.text())
        logger.info('visualizing ' + res_text)
        waveform = self.waveform_results.get_waveform(self.last_array, self.last_signal, res_text)
        if waveform:
            self.mpl_waveform.plot_waveform(waveform)
        else:
            logger.error('could not find waveform for: ' + res_text)

if __name__ == '__main__':
    # create the app
    sl_app = QtWidgets.QApplication(sys.argv)
    main = Main('config.ini')
    main.show()
    sys.exit(sl_app.exec_())
