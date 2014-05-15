__author__ = 'phoetrymaster'

from osgeo import gdal
from osgeo.gdalconst import *
import os
import sys
import subprocess

gdal.UseExceptions()

def find_files(searchdir, ext):
    foundfiles = []

    for root, dirs, files in os.walk(searchdir):
        for f in files:
            if f.upper().endswith(ext.upper()):
                foundfile = os.path.join(root, f)
                foundfiles.append(foundfile)

    return foundfiles

def create_output_raster(outFile, cols, rows, bands, datatype, driver):
    driver = gdal.GetDriverByName(driver)
    driver.Register()

    outds = driver.Create(outFile, cols, rows, bands, datatype)

    return outds

def get_output_params(filepath):
    image = gdal.Open(filepath, GA_ReadOnly)

    if image is None:
        raise Exception("Could not open " + filepath)

    rows = image.RasterYSize
    cols = image.RasterXSize
    band = image.GetRasterBand(1)
    bandtype = band.DataType
    geotransform = image.GetGeoTransform()
    projection2 = image.GetProjection()
    driver = image.GetDriver().ShortName
    numbands = image.RasterCount

    image = ""

    return rows, cols, bandtype, geotransform, projection2, driver, numbands

def clip_raster_to_extent(inraster, outrastername, xmin, ymin, xextent, yextent):

    rows, cols, datatype, geotransform, projection, driver, bands = get_output_params(inraster)

    outds = create_output_raster(outrastername, xextent, yextent, bands, datatype, driver)

    image = gdal.Open(inraster)

    for i in range(0, bands):
        print "\tProcessing band {0} of {1}...".format(i + 1, bands)
        band = image.GetRasterBand(1)
        outband = outds.GetRasterBand(i + 1)

        nodatavalue = band.GetNoDataValue()

        print "\t\tReading band data to array..."
        data = band.ReadAsArray(xmin, ymin, xextent, yextent)

        print "\t\tWriting band data to output band..."
        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(nodatavalue)
        outband.FlushCache()

        del data, outband, band

    outds.SetGeoTransform(geotransform)
    outds.SetProjection(projection)

    image = ""
    outds = ""


def find_rasters_and_multiple_clip_to_coords(rootdir, ext, coords, outformat, outputdir=""):
    if not outputdir:
        nooutputdir = 1
    else:
        nooutputdir = 0

    print "\nFinding rasters to clip in {0}.".format(rootdir)
    rasters = find_files(rootdir, ext)
    print "  Found {0} files.".format(len(rasters))
    rastcount = 1

    for raster in rasters:
        print "\nClipping raster {0} ({1} of {2}).".format(raster, rastcount, len(rasters))
        outdir, infile = os.path.split(raster)

        if nooutputdir:
            outputdir = outdir

        filename, inext = os.path.splitext(infile)
        rastcount += 1
        count = 1

        for coord in coords:
            outname = os.path.join(outputdir, filename + "_clip" + str(count) + inext)
            subprocess.call("gdal_translate -of {0} -srcwin {1} {2} {3} {4} {5} {6}".format(outformat, coord[0], coord[1], coord[2], coord[3], raster, outname), shell=True)
            #clip_raster_to_extent(raster, outname, coord[0], coord[1], coord[2], coord[3])
            print "\n  Clipped to coord {0} of {1}\n".format(count, len(coords))
            count += 1

    print "\nProcess completed."


########## PROCEDURE ##########

if __name__ == '__main__':

    ##Set Args##
    root = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected"
    ext = ".tif"
    coordslist = ((2482, 1089, 100, 100), (2892, 179, 100, 100), (3394, 1348, 100, 100), (1847, 1068, 100, 100), (1296, 355, 100, 100), (2539, 291, 100, 100))
    outdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips"
    outformat = "ENVI"

    sys.exit(find_rasters_and_multiple_clip_to_coords(root, ext, coordslist, outformat, outputdir=outdir))