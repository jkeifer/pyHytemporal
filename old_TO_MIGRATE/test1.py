from osgeo import gdal
from osgeo.gdalconst import *
import os

__author__ = 'phoetrymaster'


outDIR = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2012/"
outName = "test2.tif"
outFile = os.path.join(outDIR, outName)

inHDF = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2012/e4ftl01.cr.usgs.gov/MODIS_Composites/MOLT/MOD13Q1.005/2011.12.19/MOD13Q1.A2011353.h09v05.005.2012005014916.hdf"

hdf = gdal.Open(inHDF, GA_ReadOnly)

if hdf is None:
    print 'Could not open ' + inHDF
    sys.exit(1)

hdfDS = hdf.GetSubDatasets()

for data in hdfDS:
    if "NDVI" in data[0]:
        ndviPath = data[0]
    elif "EVI" in data[0]:
        eviPath = data[0]

hdf = ""

driver = gdal.GetDriverByName('GTiff')
driver.Register()

ndvi = gdal.Open(ndviPath)

rows = ndvi.RasterYSize
cols = ndvi.RasterXSize

band = ndvi.GetRasterBand(1)
bandType = band.DataType

data = band.ReadAsArray(0, 0, cols, rows)

outDS = driver.Create(outFile, cols, rows, 1, bandType)
outBand = outDS.GetRasterBand(1)

data = band.ReadAsArray(0, 0, cols, rows)
outBand.WriteArray(data, 0, 0)

outBand.SetNoDataValue(-99)
outBand.FlushCache()

outDS.SetGeoTransform(ndvi.GetGeoTransform())
outDS.SetProjection(ndvi.GetProjection())