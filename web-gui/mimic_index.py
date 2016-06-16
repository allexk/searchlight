from mimic import MIMICWaveformData, Patient, PatientRecord
import os
import cPickle as pickle
from collections import OrderedDict
import datetime


class MIMICIndex(object):
    def __init__(self, mimic_cache, index_file_name=None):
        """Create the MIMIC index from the MIMIC cache.

        Tries to open the existing MIMIC index if does exist. If not, creates it.
        :param mimic_cache: path to the MIMIC cache (pickled)
        :param index_file_name: index file name; if None use default
        :return: MIMIC index object
        """
        if os.path.exists(mimic_cache):
            print 'Found cached catalog in mimic.cache...'
            with open(mimic_cache, 'rb') as cache:
                self.mimic = pickle.load(cache)
        # try to find the index
        if index_file_name is None:
            index_file_name = os.path.splitext(mimic_cache)[0] + ".index"
        if os.path.exists(index_file_name):
            with open(index_file_name, "rb") as index:
                self.index = pickle.load(index)
        else:
            # need to create it and store
            self._create_index()
            print "Storing index into", index_file_name
            with open(index_file_name, "wb") as index:
                pickle.dump(self.index, index)

    def _create_index(self):
        """Create the MIMIC index: global_id -> (patient_id, record_id)."""
        self.index = dict()
        for patient in self.mimic:
            patient_id = patient.id
            for rec_id, rec in enumerate(patient):
                global_id = rec.global_id
                self.index[global_id] = (patient_id, rec_id)

    def find_segment(self, global_id, tick):
        """Find the appropriate waveform segment for the global record id and tick."""
        (patient_id, record_id) = self.index[global_id]
        print "Found record: %s-%d" % (patient_id, record_id)
        (record_name, start_time) = self.mimic.get_segment_info(patient_id, record_id, tick)
        return patient_id, record_name, start_time.strftime('%H:%M:%S.%f')[:-3] # need milliseconds


def _main():
    """Runs if the module is executed directly as a script."""
    import argparse
    parser = argparse.ArgumentParser(description=
                                     "Take the MIMIC cache and creates/uses additional indexes to find patient records")
    parser.add_argument('mimic_cache', help='MIMIC cache file')
    parser.add_argument('global_id', type=int, help='Global record id')
    parser.add_argument('tick', type=int, help='Record tick')
    opts = parser.parse_args()

    index = MIMICIndex(opts.mimic_cache)
    print "Start time is:", index.find_segment(opts.global_id, opts.tick)

if __name__ == '__main__':
    import sys
    sys.exit(_main())
