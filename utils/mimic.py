#!/usr/bin/python

import os
import sys
import re
import datetime
import scipy.io as sio
import csv
from collections import OrderedDict

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
                    segment['Id'] = i # some segments might be empty, so we need to id them
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

    def store_to_csv(self, record_mat_dir, csv_writer_data, csv_writer_segments, req_signals):
        """Store all segments signals to the specified CSV writer.

        The assumed schema is: global_rec_id, tick, <signals>. <signals> means all (currently 21) signals
         where the missing ones from the segment are marked as NULLs.

        Meta information about segments consists of the global id of the record, signal number,
         start time (in datetime and ticks) and length in ticks.

        Params:
            record_mat_dir -- directory with the MATLAB data files for the segments
            csv_writer_data -- CSV writer object capable of writing multiple rows with 'writerows()' (for data)
            csv_writer_segments -- the same as above, but for meta info
            req_signals -- signals to process
        """
        for seg in self._segments:
            csv_writer_segments.writerow([self.global_id, seg['Id'], seg['StartTime'].strftime('%H:%M:%S.%f'),
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
                        signals[tuple_idx] = '{:.3f}'.format(signals_row[array_idx])
                    csv_row = [self.global_id, start_tick] + signals
                    start_tick += 1
                    csv_writer_data.writerow(csv_row)
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
        return "Record %s with %d segments" % (self.name, len(self._segments))

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
    def __init__(self, patient_id, records, orig_patient_dir):
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
            self._records.append(PatientRecord(group_name, len(self._records), group_ord,
                                               os.path.join(orig_patient_dir, self.id)))

    def store_to_csv(self, patient_mat_dir, csv_writer_data, csv_writer_records, csv_writer_segments, req_signals):
        """Store all records to the specified CSV writer.

        For the schema see the PatientRecord's method with the same name. All records are stored one after
         another, without gaps in the CSV,

        The meta about records consists of the global record id, patient id and record name for each record.

        The meta about segments is described in the corresponding function in PatientRecord.

        Params:
            patient_mat_dir -- directory with the MATLAB data files for the patient
            csv_writer_data -- CSV writer object capable of writing multiple rows with 'writerows()'
            csv_writer_records -- the same for records meta
            csv_writer_segments -- the same for segments meta
            req_signals -- signals to process
        """
        for rec in self._records:
            csv_writer_records.writerow([rec.global_id, self.id, rec.name])
            rec.store_to_csv(patient_mat_dir, csv_writer_data, csv_writer_segments, req_signals)

    def __iter__(self):
        """Iterates over patient records."""
        for rec in self._records:
            yield rec


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
            self._patients[patient_id] = Patient(patient_id, record_files, mimic_orig_dir)

    def store_to_csv(self, output_path, req_signals):
        """Store all segments for all patients into a CSV file.

        For the schema see the PatientRecord's method with the same name. All records are stored one after
         another, without gaps in the CSV,

        Params:
            output_path -- file path for the CSV file (will be created/overwritten)
            req_signals -- signals to process
        """
        with open(os.path.join(output_path, 'mimic_wave.csv'), 'w') as csv_data_file:
            csv_writer_data = csv.writer(csv_data_file)
            csv_writer_data.writerow(['RecordId', 'Tick'] + req_signals)
            with open(os.path.join(output_path, 'mimic_records.csv'), 'w') as csv_records_file:
                csv_writer_records = csv.writer(csv_records_file)
                csv_writer_records.writerow(['RecordId', 'PatientId', 'Record'])
                with open(os.path.join(output_path, 'mimic_segments.csv'), 'w') as csv_segments_file:
                    csv_writer_segments = csv.writer(csv_segments_file)
                    csv_writer_segments.writerow(['RecordId', 'SegmentId', 'StartTime', 'StartTick', 'Ticks'])

                    for (patient_id, patient) in self._patients.iteritems():
                        patient_mat_dir = os.path.join(self._mimic_mat_dir, patient_id)
                        patient.store_to_csv(patient_mat_dir, csv_writer_data, csv_writer_records,
                                             csv_writer_segments, _list_to_ordered_dict(req_signals))

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
    parser.add_argument('--csv', metavar='path to CSV', default=None, help="Path to store signals/meta CSV")
    parser.add_argument('--signals', metavar='signals', nargs='*', default=_SIGNALS, help="Signals to process")
    parser.add_argument('mimic_mat_dir', help='Directory with MIMIC waveform MATLAB files')
    parser.add_argument('mimic_orig_dir', help='Directory with MIMIC original data/header files')

    # parse
    opts = parser.parse_args()

    # Create the catalog
    try:
        mimic = MIMICWaveformData(opts.mimic_mat_dir, opts.mimic_orig_dir)
    except ValueError:
        print "Error when parsing MIMIC directory tree:", sys.exc_info()
        return 1

    # iterate over groups of records
    all_signals = set()
    for patient  in mimic:
        print "Found following groups for patient %s:" % patient.id
        for rec in patient:
            print str(rec)
            print "Segments:"
            rec.print_segments()
            for seg in rec:
                all_signals.update(seg['Signals'])

    # check if we're missing some in the _SIGNALS
    missing_signals = set(all_signals) - set(_SIGNALS)
    if missing_signals:
        print 'Signals missing from _SIGNALS:', ', '.join(missing_signals)
    else:
        print 'No signals are missing from _SIGNALS'

    # check if we need to save the CSV
    if opts.csv:
        mimic.store_to_csv(opts.csv, opts.signals)
    return 0

if __name__ == '__main__':
    sys.exit(_main())
