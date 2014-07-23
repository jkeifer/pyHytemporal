import os
import numpy
from osgeo import gdal
from osgeo.gdalconst import *
from pyhytemporal.utils import band_number_to_doy
from pyhytemporal.vectorFunctions import get_px_coords_from_points
from utils import find_files
from core import signatureCollection


def get_sigs_in_dir(directory, viname=None, searchstring='mean.ref', recursivesearch=False):
    signatures = signatureCollection(viName=viname)
    sigFiles = find_files(directory, searchstring, recursive=recursivesearch)

    for f in sigFiles:
        signatures.add(f)

    return signatures


def get_crop_pixel_values(imagepath, locations):

    #TODO docstring

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
        for i in range(bands):
            band = img.GetRasterBand(i + 1)
            v = int(band.ReadAsArray(location[1], location[0], 1, 1))
            print "\tBand {0}: {1}".format(i + 1, v)
            values.append(v)
            band = None
        refs.append(values)
    img = None

    return refs


def write_refs_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing pixel curves to output file:"
    print referencevalues

    #TODO: Rewrite this and next function with append on open, and integrate together to reduce redundant code

    #output individual pixel curves to file
    with open(os.path.join(outdir, cropname + postfix + "_points.ref"), "w") as f:
        if comment:
            f.write("//" + comment + "\n\n")
        point = 1
        for points in referencevalues:
            f.write("\nPoint {0}:\n".format(point))
            point += 1
            imgnumber = 1
            for val in points:
                doy = band_number_to_doy(imgnumber, startdoy, doyinterval)
                f.write("{0} {1}\n".format(doy, val))
                imgnumber += 1


def write_mean_ref_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):

    #TODO docstring

    print "Writing mean reference curve to output file:"
    print referencevalues

    #output mean pixel values to file
    with open(os.path.join(outdir, cropname + postfix + "_mean.ref"), "w") as f:
        if comment:
            f.write("//" + comment + "\n\n")
        meanvals = get_mean_values(referencevalues)
        imgnumber = 1
        for val in meanvals:
            doy = band_number_to_doy(imgnumber, startdoy, doyinterval)
            f.write("{0} {1}\n".format(doy, val))
            imgnumber += 1


def get_mean_values(referencevalues):

    #TODO docstring

    array = numpy.array(referencevalues)
    mean = numpy.mean(array, axis=0)
    return list(mean)


def get_reference_curves(image, refstoget, startdoy, imageinterval, outdir="", filepostfix=""):
    """

    """

    #TODO docstring

    if not outdir:
        outdir = os.path.dirname(image)

    for shapefile in refstoget:
        cropname = os.path.splitext(os.path.basename(shapefile))[0]
        locs = get_px_coords_from_points(image, shapefile)
        refvals = get_crop_pixel_values(image, locs)
        comment = "Generated from {0} by get_crop_pixel_values version 0.1.".format(image)
        write_refs_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)
        write_mean_ref_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment,
                              postfix=filepostfix)


