__author__ = 'phoetrymaster'

from osgeo import gdal
from osgeo.gdalconst import *
import sys

def get_reference_values(imagepath, refstoget):

    from math import floor

    gdal.AllRegister()
    img = gdal.Open(imagepath, GA_ReadOnly)
    if img is None:
        raise Exception("Could not open " + imagepath)
    bands = img.RasterCount
    print "Found {0} bands in input image.".format(bands)
    refs = {}
    for key, val in refstoget.items():
        print "Processing {0} coordinates:".format(key)
        dict = {}
        for i in range(0, bands):
            band = img.GetRasterBand(i + 1)
            print "\tProcessing band {0}".format(i + 1)
            values = []
            for loc in val:
                print "\t\tGetting position {0}".format(loc)
                values.append(int(band.ReadAsArray(int(floor(loc[0])), int(floor(loc[1])), 1, 1)))
            dict[(i * 16 + 1)] = sum(values) / float(len(values))
            band = None
        refs[key] = dict
    img = None
    comment = "Generated from {0} by get_ref_values version 0.1.".format(imagepath)
    return refs, comment


def get_reference_values(imagepath, locations):

    from math import floor

    gdal.AllRegister()
    img = gdal.Open(imagepath, GA_ReadOnly)

    if img is None:
        raise Exception("Could not open " + imagepath)

    bands = img.RasterCount
    print "Found {0} bands in input image.".format(bands)

    for location in locations:
        print "Processing coordinates: {0}".format(location)
        values = []
        for i in range(0, bands):
            band = img.GetRasterBand(i + 1)
            v = (int(band.ReadAsArray(int(floor(location[0])), int(floor(location[1])), 1, 1)))
            print "\tBand {0}: {1}".format(i + 1, v)
            values.append(v)
            band = None
    img = None

    comment = "Generated from {0} by get_ref_values version 0.1.".format(imagepath)
    return refs, comment


def get_sort_dates_values(vals, threshold=-3000):
    """Gets the DOY dates (the keys) in a list from dictionary values and sorts those, placing them in chronological order
    (list x0). Then the function iterates over these values and gets the corresponding values, thresholding values if
    they are lower than an optional threshold value (-3000 default == nodata in MODIS imagery), then appending them to
    the list y. x and y are then returned."""

    x = vals.keys()
    x.sort()
    y = []

    for i in x:
        if vals[i] < threshold:
            y.append(threshold)
        else:
            y.append(vals[i])

    return x, y


def write_refs_to_txt(refs, outdir, comment=""):
    print "Writing reference curves to output files"
    for ref in refs:
        with open(os.path.join(outdir, ref[0] + ".ref"), "w") as f:
            if comment:
                f.write("//"+comment+"\n")
            keys, values = get_sort_date_values(ref[1])
            for key, val in keys, values:
                f.write("{0} {1}".format(key, val))



if __name__ == '__main__':

    image = "/Volumes/J_KEIFER/Thesis/Data/ARC_Testing/test1.dat"

    soylocs = [(6002, 2143), (5944, 2102), (5746, 2183), (5998, 2171)]
    cornlocs = [(5997, 2139), (5940, 2096), (6051, 2230), (5691, 1998)]
    wheatlocs = [(5993, 2136), (5937, 2080), (5935, 2076), (5921, 2217)]
    refstoget = {"soy": soylocs, "corn": cornlocs, "wheat": wheatlocs}

    sys.exit(get_reference_values(image, refstoget))

