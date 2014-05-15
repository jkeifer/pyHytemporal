from osgeo import gdal
from osgeo.gdalconst import *
import os, sys

__author__ = 'phoetrymaster'

rootDIR = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2012/"
#rootDIR = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS 7_2012-2013/"
outName = "test"

newfoldername = "kansas"

find = "EVI"
ext = ".hdf"
drivercode = "ENVI"
ndvalue = -3000
projection = "PROJCS[\"Sinusoidal\",GEOGCS[\"GCS_Undefined\",DATUM[\"D_Undefined\",SPHEROID[\"User_Defined_Spheroid\",6371007.181,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"False_Easting\",0.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",0.0],UNIT[\"Meter\",1.0]]"


########## METHODS ##########


def find_files(searchdir, ext):
    hdfs = []

    for root, dirs, files in os.walk(searchdir):
        for f in files:
            if f.upper().endswith(ext.upper()):
                foundfile = os.path.join(root, f)
                hdfs.append(foundfile)

    return hdfs


def create_output_dir(root, name):
    dirpath = os.path.join(root, name)

    if os.path.isdir(dirpath):
        count = 1
        dirpath_ = dirpath + "_"
        while 1:
            dirpath = dirpath_ + str(count)
            count += 1
            if not os.path.isdir(dirpath):
                break

    os.makedirs(dirpath)
    return dirpath


def create_output_raster(outFile, cols, rows, bands, datatype, drivername="GTiff"):
    driver = gdal.GetDriverByName(drivername)
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

    image = ""

    return rows, cols, bandtype, geotransform, projection


def get_hdf_subdatasets(hdfpath):
    hdf = gdal.Open(hdfpath, GA_ReadOnly)

    if hdf is None:
        raise Exception("Could not open " + hdfpath)

    sds = []
    hdfsds = hdf.GetSubDatasets()

    for data in hdfsds:
        sds.append((data[0], data[0].split(" ")[-1]))

    hdf = ""

    return sds


def main():
    outdir = create_output_dir(rootDIR, newfoldername)
    print "\nOutputting files to : {0}".format(outdir)

    print "\nFinding HDF files in directory/subfolders: {0}".format(rootDIR)
    hdfs = find_files(rootDIR, ext)
    print "\tFound {0} files.".format(len(hdfs))

    print "\nGetting images to process of type {0}...".format(find)
    toprocess = []

    for hdf in hdfs:
        sds = get_hdf_subdatasets(hdf)
        for ds in sds:
            if find.upper() in ds[1].upper():
                toprocess.append(ds[0])
                print "\t\t{0}".format(ds[0])

    bands = len(toprocess)
    print "\tFound {0} images of type {1}.".format(bands, find)

    print "\nGetting output parameters..."
    rows, cols, datatype, geotransform, projection = get_output_params(toprocess[0])
    print "\tParameters: rows: {0}, cols: {1}, datatype: {2}, projection: {3}.".format(rows, cols, datatype, projection)

    outfile = os.path.join(outdir, outName) + ".tif"
    print "\nOutput file is: {0}".format(outfile)

    outds = create_output_raster(outfile, cols, rows, bands, datatype, drivername=drivercode)
    print "\tCreated output file."

    print"\nAdding bands to output file..."
    for i in range(0, bands):
        print "\tProcessing band {0} of {1}...".format(i + 1, bands)
        image = gdal.Open(toprocess[i])
        band = image.GetRasterBand(1)

        outband = outds.GetRasterBand(i + 1)

        print "\t\tReading band data to array..."
        data = band.ReadAsArray(0, 0, cols, rows)

        print "\t\tWriting band data to output band..."
        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(ndvalue)
        outband.FlushCache()

        del data, outband
        image = ""

    print "\tFinished adding bands to output file."

    print "\nSetting transform and projection..."
    outds.SetGeoTransform(geotransform)
    outds.SetProjection(projection)

    outDS = ""

    print "\nProcess completed."


########## PROCEDURE ##########


if __name__ == '__main__':
    sys.exit(main())