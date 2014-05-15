__author__ = 'phoetrymaster'
import accuracy_assessment
import numpy as np
import os
from osgeo import gdal
from osgeo.gdalconst import *
import sys


searchdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/test1_envicurves/fullpxonly/clip1refs/KansasNDVI_2012_clip1_SLSQP/"
outFile = os.path.join(searchdir, "classified_3crop.tif")
cropimgpath = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
searchstringsvals = [("soy.", 5), ("wwheat.", 24), ("corn", 1)]#, ("sorghum", 4), ("wwheatsoydbl", 26)]
nodata = -3000

accuracyreport = os.path.join(searchdir, "accuracy_3crop.txt")

try:
    gdal.AllRegister()
    #np.set_printoptions(threshold=np.nan)

    #Crop image is constant for all iterations
    cropimg = gdal.Open(cropimgpath, GA_ReadOnly)
    if cropimg is None:
        raise Exception("Could not open: {0}".format(cropimgpath))
    else:
        rows = cropimg.RasterYSize
        cols = cropimg.RasterXSize
        projection = cropimg.GetProjection()
        transformation = cropimg.GetGeoTransform()
        band = cropimg.GetRasterBand(1)
        datatype = band.DataType
        croparray = band.ReadAsArray(0, 0, cols, rows)
        band =""
        cropimg = ""
        print "Opened crop img"

    #
    filelist = []
    files = os.listdir(searchdir)
    for f in files:
        for string, val in searchstringsvals:
            if f.endswith(".tif"):
                if string in f:
                    filelist.append((os.path.join(searchdir, f), val))


accuracy, classification, cols, rows, outstring = classify_with_threshold(croparray,
                                                                                           filelist,
                                                                                           searchdir, searchstringsvals,
                                                                                           bestthresh, nodata)
driver = gdal.GetDriverByName("ENVI")
driver.Register()

outds = driver.Create(outFile, cols, rows, 1, GDT_Int16)
outds.SetGeoTransform(transformation)
outds.SetProjection(projection)
outband = outds.GetRasterBand(1)
outband.WriteArray(classification, 0, 0)
outband.SetNoDataValue(-3000)
outband.FlushCache()

outband = ""
outds = ""

print "outputted"