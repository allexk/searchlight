#!/usr/bin/python
from __future__ import print_function

import os
import argparse
from scidb_http import SciDBConnection

# The insert command
INSERT_AFL = "insert(redimension(input(" \
             "<record: int64, tick: int64, ABP: double null, MCL1: double null, PAP: double null>" \
             "[i=0:*,10000,0], %upload%, 0, '(int64, int64, double null, double null, double null)')," \
             "mimic), mimic)"

# parse command line parameters
parser = argparse.ArgumentParser(
    description="Loads binary MIMIC data into SciDb using the HTTP interface",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--port', metavar='N', type=int, default=8080, help='Port to connect to')
parser.add_argument('--user', metavar='name', default='SciDBUser', help='User to connect as')
parser.add_argument('--password', metavar='path', help='Path to a file with the password')
parser.add_argument('host', help='Host to connect to')
parser.add_argument('patients', help='File containing list of patients')
opts = parser.parse_args()

if not opts.password:
    raise ValueError('--password is mandatory for the http method')
with open(opts.password, 'r') as passw_file:
    password = passw_file.read().strip()
scidb = SciDBConnection(opts.host, opts.port, opts.user, password)

try:
    scidb.connect()
    with open(opts.patients, 'r') as patients_file:
        for patient in patients_file:
            patient = patient.strip()
            if not patient:
                continue
            print("Loading waveform for patient ", patient)
            # load the corresponding data
            wave_file = os.path.join(patient, 'mimic_wave.dat')
            wave_file_size = os.stat(wave_file).st_size
            if not wave_file_size:
                print('File ', wave_file, ' is empty, skipping...')
                continue
            remote_file_name = scidb.upload(wave_file)
            # then, insert the data
            query = INSERT_AFL.replace('%upload%', "'%s'" % remote_file_name)
            scidb.query(query, False)
finally:
    scidb.close()
