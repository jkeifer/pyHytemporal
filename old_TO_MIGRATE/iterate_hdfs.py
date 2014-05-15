import build_multiband_image as b
from datetime import datetime as dt
import multiprocessing

start = dt.now()
print "Starting at {0}.".format(start)

maxprocesses = 1

#k2006 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2006"
#k2007 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2007"
#k2008 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2008"
k2009 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2009"
k2010 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2010"
k2011 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2011"
k2012 = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/2012"

#dict = {"2006": k2006, "2007": k2007, "2008": k2008, "2009": k2009, "2010": k2010, "2011": k2011, "2012": k2012}
dict = {"2009": k2009, "2010": k2010, "2011": k2011, "2012": k2012}

fileprefix = "Kansas"
#VIs = ["EVI", "NDVI"]
VIs = ["NDVI"]

driver = "ENVI"

nodata = -3000

processes = []
queue = []

for VI in VIs:
    print "\nProcessing {0}".format(VI)
    for key, val in dict.items():
        print "  Year: {0}\n".format(key)
        queue.append((val, fileprefix + VI + "_" + key, fileprefix + VI + "_" + key, VI, driver, nodata))

processcount = 0

for item in queue:
    p = multiprocessing.Process(target=b.build_multiband_image, args=item)
    p.start()
    processes.append(p)
    processcount += 1

    if processcount == maxprocesses:
        for p in processes:
            p.join()
        processcount = 0
        processes = []


print "PROCESS FINISHED."
print "Time elapsed: {0}".format(dt.now() - start)