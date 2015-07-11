#!/usr/bin/python

import subprocess
import shlex
import glob
import pyfits
import sys
import os

"""
This script retrieves SDSS photometry data from the specified SDSS server.
Since SDSS consists of a large number of FITS files, the script wgets the data
run by run. After each run, the required information is fetched from FITS files
and then the files are deleted. The information is output in CSV files. The
original folder structure is preserved:
    <run>/<camcol>/*.fits --> <run>/<camcol>/data.csv

For now we retrieve: ra, dec, objid, psf ugriz, model ugriz.

Caution: the sccript should be used with care. SDSS photometric data is about
3.4TB. If new attributes are required, the catalog has to be retrieved and
processed again, which puts additional burden on the SDSS server. For large
data transfers they suggest to contact the help desc.
"""

# PhotoObk base URL (re-run 301)
SDSS_BASE_URL='http://mirror.sdss3.org/sas/dr10/boss/photoObj/301'

def wget_run(target_dir, run):
    """
    Download FITS files for the specified run to the target directory.
    
    The data is retireved for all camcols (1-6). The directory structure
    under targe_dir is: <run>/<camcol>/*.fits.
    """
    wget_base_cmd = 'wget -nv -r -np -nH --reject="index.html*" --cut-dirs=5 %s/%d/%d/'
    for camcol in range(1, 7):
        wget_cmd = wget_base_cmd % (SDSS_BASE_URL, run, camcol)
        args = shlex.split(wget_cmd)
        subprocess.check_call(args, cwd=target_dir)

def pyfits_row_to_str(row):
    """
    Convert FITS table row to string. First, the required data is retrieved
    from the row and only then the string is returned.
    
    For a None row it returns the header.
    
    NOTE: this is a good place to specify what data is required
    """
    if not row:
        return 'ra,dec,objid,psfMag_u,psfMag_g,psfMag_r,psfMag_i,psfMag_z,'\
            'modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z'
    else:
        ra = str(row[72])
        dec = str(row[73])
        objid = str(row[0])
        psf_ugriz = [str(x) for x in row[94]]  # take from a FITS array
        model_ugriz = [str(x) for x in row[110]]  # take from a FITS array
    
        attrs_list = [ra, dec, objid] + psf_ugriz + model_ugriz
        return ','.join(attrs_list)  # join CSV style


def process_run_dir(run_path):
    """
    Process run directory and create CSV with the required information.
    
    Each camcol is processed and one CSV per camcol is created.
    """
    for camcol in range(1, 7):
        camcol_dir = os.path.join(run_path, str(camcol))
        print 'Processing directory ' + camcol_dir
        fits_files = glob.glob(os.path.join(camcol_dir, '*.fits'))
        if not fits_files:
            continue

        with open(os.path.join(camcol_dir, 'data.csv'), 'w') as out_file:
            # header
            out_file.write(pyfits_row_to_str(None))
            out_file.write('\n')
            # data
            for ff in fits_files:
                f = pyfits.open(ff)
                # SDSS photometry contains two section: the table is the second
                try:
                    data = f[1].data
                    for row in data:
                        row_str = pyfits_row_to_str(row)
                        out_file.write(row_str)
                        out_file.write('\n')
                except:
                    # some FITS don't have binary tables...
                    # we just ignore all fancy FITS data
                    pass

                # we also delete the file (don't have enough space...)
                f.close()
                os.unlink(ff)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: sdss_wget.py <spec_file>"
        sys.exit(1)

    target_dir = os.getcwd()
    with open(sys.argv[1], 'r') as spec_file:
        for run in spec_file:
            run = run.strip()
            if not run:
                continue

            run_path = os.path.join(target_dir, run)
            if not os.path.exists(run_path):
                print 'Retrieveing remote data for ' + run_path
                wget_run(target_dir, int(run))
            print 'Processing run data in ' + run_path
            process_run_dir(run_path)

    sys.exit(0)
