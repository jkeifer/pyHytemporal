__author__ = 'phoetrymaster'

import subprocess
import os
from datetime import datetime as dt
import multiprocessing

def find_files(searchdir, ext):
    foundfiles = []

    for root, dirs, files in os.walk(searchdir):
        for f in files:
            if f.upper().endswith(ext.upper()):
                foundfile = os.path.join(root, f)
                foundfiles.append(foundfile)

    return foundfiles


def reproject_modis(pathtoresampletool, inhdf, outfile, parameters):
    subprocess.call("{0} -p \"{1}\" -i \"{2}\" -o \"{3}\"".format(pathtoresampletool, parameters, inhdf, outfile), shell=True)




if __name__ == '__main__':

    start = dt.now()
    print start

    pathtoresampletool = "/Users/phoetrymaster/Applications/ModisTools/MRT/bin/resample"
    outfolder = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/"
    parameterfile = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/projectparamters.prm"
    searchdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012"

    workers = 3

    print "\nSearching for .hdf files in {0}".format(searchdir)
    hdfs = find_files(searchdir, ".hdf")
    print "\nFound {0} .hdf files to convert.\n".format(len(hdfs))

    count = 1
    processes = []
    for hdf in hdfs:
        print "Reprojecting file {0} of {1}: {2}.".format(count, len(hdfs), hdf)
        p = multiprocessing.Process(target=reproject_modis, args=(pathtoresampletool, hdf, os.path.join(outfolder, os.path.basename(hdf)), parameterfile))
        p.start()
        processes.append(p)
        count += 1

        if len(processes) == workers:
            for p in processes:
                p.join()
            processes = []

    print "\n\nProcess completed."
    print "Time elapsed: {0}".format(dt.now() - start)