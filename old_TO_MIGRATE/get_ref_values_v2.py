__author__ = 'phoetrymaster'

from osgeo import gdal
from osgeo.gdalconst import *
import sys
import os
import numpy
from math import floor
from get_px_coords_from_point import get_px_coords_from_points


def get_crop_pixel_values(imagepath, locations):

    gdal.AllRegister()
    img = gdal.Open(imagepath, GA_ReadOnly)

    if img is None:
        raise Exception("Could not open " + imagepath)

    bands = img.RasterCount
    print "Found {0} bands in input image.".format(bands)

    refs = []
    for location in locations:
        print "Processing coordinates: {0}".format(location)
        values = []
        for i in range(0, bands):
            band = img.GetRasterBand(i + 1)
            v = int(band.ReadAsArray(int(floor(location[0])), int(floor(location[1])), 1, 1))
            print "\tBand {0}: {1}".format(i + 1, v)
            values.append(v)
            band = None
        refs.append(values)
    img = None

    return refs


def write_refs_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing pixel curves to output file:"
    print referencevalues

    #output individual pixel curves to file
    with open(os.path.join(outdir, cropname + postfix + "_points.ref"), "w") as f:
        if comment:
            f.write("//"+comment+"\n\n")
        point = 1
        for points in referencevalues:
            f.write("\nPoint {0}:\n".format(point))
            point += 1
            imgnumber = 0
            st = startdoy
            for val in points:
                doy = st + imgnumber * doyinterval
                if doy > 365:
                    doy = 366
                    st = 366
                    imgnumber = 0
                f.write("{0} {1}\n".format(doy, val))
                imgnumber += 1


def write_mean_ref_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing mean reference curve to output file:"
    print referencevalues

    #output mean pixel values to file
    with open(os.path.join(outdir, cropname + postfix + "_mean.ref"), "w") as f:
        if comment:
            f.write("//"+comment+"\n\n")
        meanvals = get_mean_values(referencevalues)
        imgnumber = 0
        st = startdoy
        for val in meanvals:
            doy = int(startdoy) + imgnumber * int(doyinterval)
            if doy > 365:
                doy = 366
                st = 366
                imgnumber = 0
            f.write("{0} {1}\n".format(doy, val))
            imgnumber += 1


def get_mean_values(referencevalues):
    array = numpy.array(referencevalues)
    mean = numpy.mean(array, axis=0)
    return list(mean)


def get_reference_curves(image, refstoget, startdoy, imageinterval, outdir="", filepostfix=""):
    if not outdir:
        outdir = os.path.dirname(image)

    for key, val in refstoget.items():
        if val:
            cropname = key
            locs = val
            refvals = get_crop_pixel_values(image, locs)
            comment = "Generated from {0} by get_crop_pixel_values version 0.1.".format(image)
            write_refs_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)
            write_mean_ref_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)


if __name__ == '__main__':

    image = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/ARC_Testing/test1.dat"

    soylocs = [(6002, 2143), (5944, 2102), (5746, 2183), (5998, 2171)]
    cornlocs = [(5997, 2139), (5940, 2096), (6051, 2230), (5691, 1998)]
    wheatlocs = [(5993, 2136), (5937, 2080), (5935, 2076), (5921, 2217)]
    refstoget = {"soy": soylocs, "corn": cornlocs, "wheat": wheatlocs}


    sys.exit(get_reference_curves(image, refstoget, 1, 16, "/Users/phoetrymaster/Desktop"))

