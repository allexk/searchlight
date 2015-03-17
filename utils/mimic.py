#!/usr/bin/python

import os
import sys
import re
import datetime
import scipy.io as sio
import csv
from collections import OrderedDict
import struct
import cPickle as pickle

# Signals that can be found in MIMIC
_SIGNALS = ['I', 'II', 'III', 'V', 'PAP', 'MCL1', 'ABP', 'AVF', 'MCL', 'ART', 'AOBP', 'UAP', 'RESP', 'L', 'CVP', 'AVL',
            'R', 'AVR', 'RAP', 'PLETH', 'LAP']


def _list_to_ordered_dict(list_):
    """Create an ordered dictionary for the list.

    The dictionary contains key --> list position mapping.
    """
    return OrderedDict([(s, i) for i, s in enumerate(list_)])


def _parse_file_name(file_name):
    """Parses the file name into record name and ordinal.
    :param file_name: record file name (XXXXXXX_YYYY.mat, where XXXXXXX is the name; YYYY -- ordinal)
    :return: (name, ordinal) pair
    """
    rec_parts = re.split('_|\.', file_name)
    if len(rec_parts) != 3:
        raise ValueError('Unknown record format: %s' % file_name)
    return rec_parts[0], int(rec_parts[1])


def _parse_record_name(file_name):
    """Parse out the record name from the file name.
    :param file_name: file_name (same format as _parse_file_name())
    :return: record name
    """
    return _parse_file_name(file_name)[0]


def _ticks_to_timedelta(ticks):
    """Convert time interval in ticks into proper datetime.timedelta.

    Assume each tick corresponds to 0.008 seconds, since most (all?) MIMIC signals are 125Hz.

    :param ticks: length of time interval in ticks
    :return: time interval as datetime.timedelta
    """
    secs = ticks / 125 # Standard signal is 125Hz
    milli_secs = (ticks - 125 * secs) * 8 # 125Hz = 0.008s = 8milliseconds per tick
    return datetime.timedelta(seconds=secs, milliseconds=milli_secs)


class BinarySciDBWriter(object):
    """Class for dumping data into SciDB binary format.

    For the format description see the documentation. This class uses Python struct module to
    produce correct binary strings.
    """
    TYPE_STRUCT_MAP = {'int64': 'q', 'double': 'd', 'string': 's'}

    def __init__(self, file_, row_format):
        """Create a new binary writer according to the format.

        The format must be a sequence
         conforming to the following specification: <type> [null]...

        :param file_: file object to write the data
        :param row_format: format for the rows to write
        """
        self._file = file_
        self._types = []
        for field_type in row_format:
            field_type = field_type.strip().split()
            assert field_type, 'Field type must be of format: type [null]'
            if field_type[0] not in self.TYPE_STRUCT_MAP.keys():
                raise ValueError('Type %s is not known' % field_type[0])
            if len(field_type) > 1:
                assert field_type[1] == 'null', 'The second field of the type must be null'
                self._types.append((field_type[0], True))
            else:
                self._types.append((field_type[0], False))

        # if we don't have string types, we can pre-compute the structure
        self._struct = None
        self._format_string = ['=']
        self._string_params = {}
        for i, (type_, is_null) in enumerate(self._types):
            if is_null:
                self._format_string.append('b') # null indicator
            if type_ == 'string':
                self._format_string.append(None) # will figure out dynamically
                self._string_params[i] = len(self._format_string) - 1
            else:
                self._format_string.append(self.TYPE_STRUCT_MAP[type_])

        if 'string' not in [x[0] for x in self._types]:
            self._struct = struct.Struct(''.join(self._format_string))

    def writerow(self, row):
        """Write the specified row in the binary SciDB format.
        :param row: row to write (a sequence of values)
        """
        value_seq = []
        assert len(self._types) == len(row), 'Types/Values mismatch in binary writer'
        for (i, val) in enumerate(row):
            type_ = self._types[i]
            if type_[1]:
                # must write NULL indicator
                if val is None:
                    value_seq.append(0) # always use 0 as the NULL reason code
                else:
                    value_seq.append(-1)
            if type_[0] == 'string':
                string_format = 'i'
                if val is None:
                    assert type_[1], 'NULL string is given for the non-NULL type'
                    value_seq.append(0)
                else:
                    str_len = len(val) + 1  # always 0-terminated as requested by SciDB
                    value_seq.append(str_len)
                    string_format += '%ds' % str_len
                    value_seq.append(val)  # struct will 0-terminate it due to padding
                self._format_string[self._string_params[i]] = string_format
            else:
                # non-string types
                if val is None:
                    val = 0 # padding for NULL values (required by SciDB)
                value_seq.append(val)
        # serialize
        if self._struct:
            buf = self._struct.pack(*value_seq)
        else:
            buf = struct.pack(''.join(self._format_string), *value_seq)
        # write
        self._file.write(buf)


class PatientRecord(object):
    """Represents a record for a single patient.

    A record might consist of several files, each with the same name, but different ids.
    """

    def __init__(self, record_name, global_record_id, files_num, orig_patient_dir):
        """Creates a new record for a patient.

        A record consists of a name and the number of files in the record.
        """
        self.name = record_name
        self.global_id = global_record_id

        print "Parsing meta for record:", self.name
        self._segments = []
        tick_count = 0
        possible_formats = ['%H:%M:%S.%f', '%M:%S.%f', '%H:%M:%S']
        for i in range(1, files_num + 1):
            header_name = os.path.join(orig_patient_dir, '%s_%04d.hea' % (self.name, i))
            try:
                with open(header_name, 'r') as header:
                    segment = {}
                    meta = header.readline().strip().split()
                    # Time might come in different formats. Try them all in succession...
                    for format_ in possible_formats:
                        try:
                            segment['StartTime'] = datetime.datetime.strptime(meta[-1], format_)
                            break
                        except ValueError:
                            pass
                    if 'StartTime' not in segment:
                        raise ValueError("Cannot parse time from header:"), header_name

                    segment['Ticks'] = int(meta[-2])
                    segment['StartTick'] = tick_count
                    segment['Id'] = i  # some segments might be empty, so we need to id them
                    tick_count += segment['Ticks']

                    # A file might contain several signals
                    segment['Signals'] = []
                    signals_num = int(meta[1])
                    for s in range(signals_num):
                        signal_meta = header.readline().strip().split()
                        segment['Signals'].append(signal_meta[-1])

                    # New segment
                    self._segments.append(segment)
            except IOError:
                # No header? Just consider the corresponding segment empty.
                print "Cannot find %s, assume the segment is empty..." % header_name

    def record_length(self):
        """Return the total length of the record."""
        if not self._segments:
            return 0
        else:
            last_segment = self._segments[-1]
            return last_segment['StartTick'] + last_segment['Ticks']

    def store_to_file(self, record_mat_dir, writer_data, writer_segments, req_signals):
        """Store all segments signals to the specified writer.

        The assumed schema is: global_rec_id, tick, <signals>. <signals> means all (currently 21) signals
         where the missing ones from the segment are marked as NULLs.

        Meta information about segments consists of the global id of the record, signal number,
         start time (in datetime and ticks) and length in ticks.

        Params:
            record_mat_dir -- directory with the MATLAB data files for the segments
            writer_data -- writer object capable of writing multiple rows with 'writerows()' (for data)
            writer_segments -- the same as above, but for meta info
            req_signals -- signals to process
        """
        for seg in self._segments:
            writer_segments.writerow([self.global_id, seg['Id'], seg['StartTime'].strftime('%H:%M:%S.%f'),
                                      seg['StartTick'], seg['Ticks']])
            signal_indexes = []
            for (i, sig) in enumerate(seg['Signals']):
                assert sig in _SIGNALS, 'Unknown signal found: %s' % sig
                if sig in req_signals:
                    signal_indexes.append((req_signals[sig], i))
            if not signal_indexes:
                # this segment doesn't contain the required signals
                continue
            seg_file_name = os.path.join(record_mat_dir, '%s_%04d.mat' % (self.name, seg['Id']))
            try:
                segment_data = sio.loadmat(seg_file_name)
                # Iterate over rows of the signal array
                start_tick = seg['StartTick']
                for signals_row in segment_data['signal']:
                    signals = [None] * len(req_signals.keys())
                    for tuple_idx, array_idx in signal_indexes:
                        signals[tuple_idx] = signals_row[array_idx]
                    data_row = [self.global_id, start_tick] + signals
                    start_tick += 1
                    writer_data.writerow(data_row)
            except IOError:
                    print 'Cannot open %s to read signals, assuming the segment empty...' % seg_file_name

    def print_segments(self):
        for seg in self._segments:
            len_ticks = seg['Ticks']
            print '\t[%d]: Start: %s, Length (ticks/time): %d/%s, StartTick: %d, Signals: %s' % (
                seg['Id'], seg['StartTime'].strftime('%H:%M:%S.%f'), len_ticks, _ticks_to_timedelta(len_ticks),
                seg['StartTick'], seg['Signals']
            )

    def __str__(self):
        """Create string representation."""
        return "Record %s[%d] with %d segments" % (self.name, self.global_id, len(self._segments))

    def __iter__(self):
        """Iterate over the record's segments."""
        for seg in self._segments:
            yield seg


class Patient(object):
    """Contains all records for a specific patient.

    Each patient has groups of records, which correspond to different observations.
    Records in each group are supposed to be more or less bundled together, time-wise.

    Users can iterate over groups via the groups() generator.
    """
    def __init__(self, patient_id, records, orig_patient_dir, first_record_id):
        self.id = patient_id
        self._records = []
        # parse records into groups
        assert records, 'Records must not be empty.'
        records.sort()
        curr_index = 0
        while curr_index < len(records):
            # start another group
            group_ord = 0
            group_name = _parse_record_name(records[curr_index])
            while curr_index < len(records):
                rec = records[curr_index]
                rec_name, rec_ord = _parse_file_name(rec)
                if rec_name == group_name:
                    group_ord += 1
                    if group_ord != rec_ord:
                        # Apparently some '.mat' files are missing
                        print 'Expected file: %s_%04d.mat, found: %s_%04d.mat' % (
                            rec_name, group_ord, rec_name, rec_ord
                        )
                        group_ord = rec_ord
                    curr_index += 1
                else:
                    break
            # Create another record
            self._records.append(PatientRecord(group_name, first_record_id + len(self._records), group_ord,
                                               os.path.join(orig_patient_dir, self.id)))

    def store_to_file(self, patient_mat_dir, writer_data, writer_records, writer_segments, req_signals):
        """Store all records to the specified writer.

        For the schema see the PatientRecord's method with the same name. All records are stored one after
         another, without gaps,

        The meta about records consists of the global record id, patient id and record name for each record.

        The meta about segments is described in the corresponding function in PatientRecord.

        Params:
            patient_mat_dir -- directory with the MATLAB data files for the patient
            writer_data -- writer object capable of writing multiple rows with 'writerows()'
            writer_records -- the same for records meta
            writer_segments -- the same for segments meta
            req_signals -- signals to process
        """
        for rec in self._records:
            writer_records.writerow([rec.global_id, self.id, rec.name])
            rec.store_to_file(patient_mat_dir, writer_data, writer_segments, req_signals)

    def __iter__(self):
        """Iterates over patient records."""
        for rec in self._records:
            yield rec

    def records_num(self):
        """Return the number of records for the patient."""
        return len(self._records)


class MIMICWaveformData(object):
    """Contains all patients with their records."""

    def __init__(self, mimic_mat_dir, mimic_orig_dir):
        """Create a new MIMIC catalog from the specified directory."""
        self._patients = OrderedDict()  # id -> Patient
        self._mimic_mat_dir = mimic_mat_dir
        self._mimic_orig_dir = mimic_orig_dir
        # first, go through all patients
        patient_dirs = os.listdir(self._mimic_mat_dir)
        patient_dirs.sort()
        if not len(patient_dirs):
            raise ValueError('No patients dirs found in %s' % self._mimic_mat_dir)

        # Assume a directory per patient
        print "Parsing MIMIC catalog..."
        global_record_id = 0
        for patient_id in patient_dirs:
            patient_dir = os.path.join(self._mimic_mat_dir, patient_id)
            if not os.path.isdir(patient_dir):
                print '%s is not a directory, skipping...' % patient_dir
                continue
            record_files = os.listdir(patient_dir)
            if not len(record_files):
                print 'No record files for: %s, skipping...' % patient_id
                continue
            # Create a new patient
            self._patients[patient_id] = Patient(patient_id, record_files, mimic_orig_dir, global_record_id)
            global_record_id += self._patients[patient_id]

    def store_to_file(self, output_path, req_signals, patients_filter, binary):
        """Store all segments for all patients into a file.

        For the schema see the PatientRecord's method with the same name. All records are stored one after
         another, without gaps,

        Params:
            output_path -- path for the file (will be created/overwritten)
            req_signals -- signals to process
            binary -- True, if use the SciDB binary format
        """
        if binary:
            files_ext = '.dat'
        else:
            files_ext = '.csv'
        with open(os.path.join(output_path, 'mimic_wave%s' % files_ext), 'w') as data_file:
            if binary:
                writer_data = BinarySciDBWriter(data_file, ['int64', 'int64'] + ['double null'] * len(req_signals))
            else:
                writer_data = csv.writer(data_file)
                writer_data.writerow(['RecordId', 'Tick'] + req_signals)
            with open(os.path.join(output_path, 'mimic_records%s' % files_ext), 'w') as records_file:
                if binary:
                    writer_records = BinarySciDBWriter(records_file, ['int64', 'string', 'string'])
                else:
                    writer_records = csv.writer(records_file)
                    writer_records.writerow(['RecordId', 'PatientId', 'Record'])
                with open(os.path.join(output_path, 'mimic_segments%s' % files_ext), 'w') as segments_file:
                    if binary:
                        writer_segments = BinarySciDBWriter(segments_file,
                                                            ['int64', 'int64', 'string', 'int64', 'int64'])
                    else:
                        writer_segments = csv.writer(segments_file)
                        writer_segments.writerow(['RecordId', 'SegmentId', 'StartTime', 'StartTick', 'Ticks'])
                    # we have created the CSV files. Now dump all patients (subject to the filter)
                    for (patient_id, patient) in self._patients.iteritems():
                        if not patients_filter or patient_id in patients_filter:
                            patient_mat_dir = os.path.join(self._mimic_mat_dir, patient_id)
                            patient.store_to_file(patient_mat_dir, writer_data, writer_records,
                                                  writer_segments, _list_to_ordered_dict(req_signals))

    def __iter__(self):
        """Iterates over all patients."""
        for patient in self._patients.values():
            yield patient


def _main():
    """Runs if the module is executed directly as a script."""
    # command line parser
    import argparse
    parser = argparse.ArgumentParser(description=
                                     "Scans the MIMIC waveform data and parses into CSV (modifying the layout)")

    # command line parameters
    parser.add_argument('--store', metavar='path to CSV', default=None, help="Path to store signals/meta CSV")
    parser.add_argument('--signals', metavar='signals', nargs='*', default=_SIGNALS, help="Signals to process")
    parser.add_argument('--filter', metavar='patients file', help="File containing patient ids to load")
    parser.add_argument('--binary', action='store_true', help='Write in the SciDB binary format')
    parser.add_argument('--show-catalog', action='store_true', help='Output the catalog to stdout')
    parser.add_argument('--single-binary', action='store_true', help='Store data in a single binary')
    parser.add_argument('mimic_mat_dir', help='Directory with MIMIC waveform MATLAB files')
    parser.add_argument('mimic_orig_dir', help='Directory with MIMIC original data/header files')

    # parse
    opts = parser.parse_args()

    # check if we have cache
    if os.path.isfile('mimic.cache'):
        print 'Found cached catalog in mimic.cache...'
        with open('mimic.cache', 'rb') as mimic_cache:
            mimic = pickle.load(mimic_cache)
    else:
        # Create the catalog
        print 'Creating MIMIC catalog...'
        try:
            mimic = MIMICWaveformData(opts.mimic_mat_dir, opts.mimic_orig_dir)
        except ValueError:
            print "Error when parsing MIMIC directory tree:", sys.exc_info()
            return 1
        # Create cache
        print 'Caching the catalog in mimic.cache...'
        with open('mimic.cache', 'wb') as mimic_cache:
            pickle.dump(mimic, mimic_cache)

    # iterate over groups of records
    if opts.show_catalog:
        max_record_length = 0
        all_signals = set()
        for patient in mimic:
            print "Found following groups for patient %s:" % patient.id
            for rec in patient:
                print str(rec)
                print "Segments:"
                rec.print_segments()
                for seg in rec:
                    all_signals.update(seg['Signals'])
                max_record_length = max(max_record_length, rec.record_length())
        print 'Maximum record length: %d' % max_record_length

        # check if we're missing some in the _SIGNALS
        missing_signals = set(all_signals) - set(_SIGNALS)
        if missing_signals:
            print 'Signals missing from _SIGNALS:', ', '.join(missing_signals)
        else:
            print 'No signals are missing from _SIGNALS'

    # check if we need to save the signals
    if opts.store:
        patients = None
        if opts.filter:
            with open(opts.filter, 'r') as filter_file:
                patients = set([l.strip() for l in filter_file if l.strip()])
        if patients and not opts.single_binary:
            # store each patient in its own dir
            for patient in patients:
                store_path = os.path.join(opts.store, str(patient))
                if not os.path.exists(store_path):
                    os.mkdir(store_path)
                else:
                    print "'%s' already exists. It will be rewritten..." % store_path
                print 'Dumping data for %s into %s...' % (patient, store_path)
                mimic.store_to_file(store_path, opts.signals, [patient], opts.binary)
        else:
            mimic.store_to_file(opts.store, opts.signals, patients, opts.binary)
    return 0

if __name__ == '__main__':
    sys.exit(_main())
