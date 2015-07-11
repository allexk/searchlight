#!/usr/bin/python

import csv
import os
import sys

# multiplier for ra/dec coordinates
RA_DEC_MULT=10000

def process_csv(path, fields, writer, min_max_ra, min_max_dec):
    with open(path, 'rb') as cf:
        print 'Processing file', path
        reader = csv.reader(cf)
        
        # assume the header goes first
        header = reader.next()
        try:
            indexes = [header.index(f) for f in fields]
        except ValueError, ex:
            print str(ex)
            raise
        
        try:
            for row in reader:
                # ra
                ra = int(float(row[0]) * RA_DEC_MULT)
                min_max_ra[0] = min(min_max_ra[0], ra)
                min_max_ra[1] = max(min_max_ra[1], ra)

                # dec
                dec = int(float(row[1]) * RA_DEC_MULT)
                min_max_dec[0] = min(min_max_dec[0], dec)
                min_max_dec[1] = max(min_max_dec[1], dec)

                # the rest
                vals = [row[i] for i in indexes]

                writer.writerow([ra, dec] + vals)
        except csv.Error, ex:
            print 'ex'


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: sdss_walker.py <sdss_csv_folder> <csv_file> <fields...>'
        sys.exit(1)
        
    # params
    fields = sys.argv[3:]
    min_max_ra = [1000, -1000]
    min_max_dec = [1000, -1000]

    # writer and header
    with open(sys.argv[2], 'wb') as csv_outf:
        writer = csv.writer(csv_outf, lineterminator='\n')
        writer.writerow(['ra', 'dec'] + fields)

        # go!
        for root, dirs, files in os.walk(sys.argv[1]):
            print 'Processing files in', root
            csv_files = [f for f in files if f.endswith('.csv')]
            for cf in csv_files:
                process_csv(os.path.join(root, cf), fields, writer,
                            min_max_ra, min_max_dec)

    print 'Stats: ra=%s, dec=%s' % (str(min_max_ra), str(min_max_dec))
    sys.exit(0)
