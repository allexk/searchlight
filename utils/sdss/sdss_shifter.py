#!/usr/bin/python

import csv
import sys

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'Usage: sdss_shifter.py <in_csv> <out_csv>'
        sys.exit(1)
        
    # reader
    with open(sys.argv[1], 'rb') as in_csv, open(sys.argv[2], 'wb') as out_csv:
        # header
        reader = csv.reader(in_csv)
        header = reader.next()

        writer = csv.writer(out_csv, lineterminator='\n')
        writer.writerow(header)

        # init max/min
        row = reader.next()
        min_max_ra = [int(row[0]), int(row[0])]
        min_max_dec = [int(row[1]), int(row[1])]

        for row in reader:
            ra = int(row[0])
            dec = int(row[1])
            dec += 250574
            row[1] = dec


            # check for weird rows
            # ra = int(row[0])
            # dec = int(row[1])
            # if ra < -10000000 or ra > 100000000:
            #    print 'Omitting row:', row
            # elif dec < -10000000 or dec > 100000000:
            #    print 'Omitting row:', row

            min_max_ra[0] = min(min_max_ra[0], ra)
            min_max_ra[1] = max(min_max_ra[1], ra)

            min_max_dec[0] = min(min_max_dec[0], dec)
            min_max_dec[1] = max(min_max_dec[1], dec)

            writer.writerow(row)


    print 'Stats: ra=%s, dec=%s' % (str(min_max_ra), str(min_max_dec))
    sys.exit(0)
