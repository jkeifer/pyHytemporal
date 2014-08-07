import os
import subprocess
import numpy
import sys
from osgeo import gdal
from osgeo.gdalconst import *
from osgeo.gdalconst import GA_ReadOnly
from osgeo import osr
from core import gdalObject, gdalProperties
from utils import change_geotransform, create_output_dir, find_files
from vectorFunctions import get_ref_from_shapefile, read_shapefile_to_points, get_px_coords_from_geographic_coords, get_geographic_coords_from_px_coords


def createNewImage(outfilepath, cols, rows, bands, datatype,
                       drivername=None, geotransform=None, projection=None):
        """
        Creates a new image file using specifed image properties.

        Required Argument(s):
            - outfilepath: The output path for the new file.
            - cols: Number of columns in the output image.
            - rows: Number of rows in the output image.
            - bands: Number of bands in the output image.
            - datatype: The datatype for the bands in the output image.

        Optional Argument(s):
            - drivername: This is the file type for the output specified by GDAL driver name. The default is GeoTIFF.
            - geotransform: This is the geotransform for the output image.
            - projection: This is the projection for the output image.

        Returns:
            - None
        """

        if drivername is None:
            drivername = "GTiff"

        driver = gdal.GetDriverByName(drivername)
        driver.Register()

        gdalImage = driver.Create(outfilepath, cols, rows, bands, datatype)

        if gdalImage is None:
            raise Exception("Error encountered creating output file.") #TODO convert to classed exception
        else:

            if geotransform:
                gdalImage.SetGeoTransform(geotransform)

            if projection:
                gdalImage.SetProjection(projection)

        return gdalImage


def copySchemaToNewImage(gdalPropertiesObject, outfilepath, numberofbands=None, drivername=None, datatype=None,
                         cols=None, rows=None, geotransform=None, projection=None):
        """
        Creates a new image using the properties of an existing image which has been loaded into a gdalPropertiesObject.

        Required Argument(s):
            - gdalPropertiesObject: A properties object containing all the original image's attributes.
            - outfilepath: The output path for the new file.

        Optional Argument(s):
            - numberofbands: If this is not specified, it will use the number of bands in the input image.
            - drivername: This is the file type for the output specified by GDAL driver name. The default is GeoTIFF.
            - datatype: The datatype of the pixels in the new image.
            - cols: The number of columns to use in the new image.
            - rows: The number of rows to use in the new image.
            - geotransform: The geotransform to set for the new image.
            - projection: The projection to set for the new image.

            Any omitted optional arguments, unless otherwise specified, will have a default value set by that property
            of the copied image.

        Returns:
            - newimage: gdal image object of the new image.
        """

        if drivername is None:
            drivername = "GTiff"

        if datatype is None:
            datatype = gdalPropertiesObject.datatype

        if cols is None:
            cols = gdalPropertiesObject.cols

        if rows is None:
            rows = gdalPropertiesObject.rows

        if geotransform is None:
            geotransform = gdalPropertiesObject.geotransform

        if projection is None:
            projection = gdalPropertiesObject.projection

        if not numberofbands:
            numberofbands = gdalPropertiesObject.bands

        newimage = createNewImage(outfilepath, cols, rows, numberofbands, datatype, drivername=drivername,
                                  geotransform=geotransform, projection=projection)

        return newimage


def openImage(infilepath, readonly=True):
    """
    Opens a raster file and updates the instance attributes.

    Required Argument(s):
        - infilepath: The path to the raster to be opened.

    Optional Argument(s):
        - None

    Returns:
        - None
    """

    gdal.AllRegister()
    if readonly is True:
        image = gdal.Open(infilepath, GA_ReadOnly)
    elif readonly is False:
        image = gdal.Open(infilepath)
    else:
        raise Exception("Error: the read status could not be be determined.")

    if image is None:
        raise Exception("Error encountered opening file.") #TODO convert to classed exception

    return image


def clip_raster_to_extent(inraster, outraster, xmin, ymin, xextent, yextent):
    """

    """
    #TODO Docstring

    original = openImage(inraster)
    originalproperties = gdalProperties(original)

    newgeotransform = change_geotransform(originalproperties.geotransform, xmin, ymin)

    outds = copySchemaToNewImage(originalproperties, outraster, cols=xextent, rows=yextent, geotransform=newgeotransform)

    for i in range(1, originalproperties.bands + 1):
        band = original.GetRasterBand(i)
        outband = outds.GetRasterBand(i)

        nodatavalue = band.GetNoDataValue()
        data = band.ReadAsArray(xmin, ymin, xextent, yextent)

        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(nodatavalue)
        outband.FlushCache()

        del data, outband, band

    original = ""
    outds = ""

    return outraster


def clip_and_mask_raster_with_shapefile(inraster, shapefile, outraster):
    """

    """
    # TODO Docstring

    from PIL import Image, ImageDraw

    raster = openImage(inraster)
    raster_properties = gdalProperties(raster)

    rasterwkt = raster.GetProjectionRef()

    oSRSop = osr.SpatialReference()
    oSRSop.ImportFromWkt(rasterwkt)

    shpextent, shppoints = read_shapefile_to_points(shapefile, oSRSop)

    points = []
    for p in xrange(shppoints.GetPointCount()):
        points.append(shppoints.GetPoint(p))

    pnts = numpy.array(points).transpose()

    cornerpnts = [(shpextent[0], shpextent[3]),  # top left
                  (shpextent[1], shpextent[2])]  # bottom right

    # TODO change this to use the function from vectorFunctions
    pixel, line = world2Pixel(raster_properties.geotransform, pnts[0], pnts[1])

    rasterPoly = Image.new("L", (raster.RasterXSize, raster.RasterYSize), 1)
    rasterize = ImageDraw.Draw(rasterPoly)
    listdata = [(pixel[i], line[i]) for i in xrange(len(pixel))]
    rasterize.polygon(listdata, 0)
    mask = 1 - PIL_image_to_array(rasterPoly)

    # find extent of new image by getting corner coords in image and subtract mins from maxes
    pxcoords = get_px_coords_from_geographic_coords(raster_properties, cornerpnts)

    ymin, xmin = pxcoords[0]
    ymax, xmax = pxcoords[1]
    xextent = xmax - xmin
    yextent = ymax - ymin


    geotransform_coords = get_geographic_coords_from_px_coords(raster_properties, [pxcoords[0]])[0]

    newgeotransform = list(raster_properties.geotransform)
    newgeotransform[0], newgeotransform[3] = geotransform_coords[0], geotransform_coords[1]

    outds = copySchemaToNewImage(raster_properties, outraster, cols=xextent, rows=yextent, geotransform=newgeotransform)

    for i in range(1, raster_properties.bands + 1):
        band = raster.GetRasterBand(i)
        outband = outds.GetRasterBand(i)

        nodatavalue = band.GetNoDataValue()
        data = band.ReadAsArray(xmin, ymin, xextent, yextent)
        data[mask[ymin:ymax, xmin:xmax] == 0] = raster_properties.nodata

        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(nodatavalue)
        outband.FlushCache()

        del data, outband, band

    raster = ""
    outds = ""

    return outraster


def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the pixel location of a geospatial coordinate
    """

    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = numpy.round((x - ulX) / xDist).astype(numpy.int)
    line = numpy.round((ulY - y) / xDist).astype(numpy.int)

    return pixel, line


def PIL_image_to_array(image):
    """
    Converts a Python Imaging Library array to a numpy array.
    """
    a = numpy.fromstring(image.tostring(),'b')
    a.shape = image.im.size[1], image.im.size[0]

    return a


def read_image_into_array(gdalimage):
    """

    """
    #TODO Docstring

    properties = gdalProperties(gdalimage)

    array = numpy.empty([properties.bands, properties.rows, properties.cols], dtype=int)

    for i in range(properties.bands):
        band = gdalimage.GetRasterBand(i + 1)
        array[i] = band.ReadAsArray(0, 0, properties.cols, properties.rows)

    array = array.transpose(1, 2, 0)  # Turn the array to allow a single pixel to be isolated in the stack (x, y, time orientation)

    return array.squeeze()  # Squeeze removes z-dimension if image only has one band


def create_test_image(imagepath, imagename, drivercode="ENVI"):
    """
    Outputs a 100px by 100px image with constant pixel values of 1000.

    "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/TestingImagery.tif"
    """

    #TODO: test for ext on imagename; add if missing -- how to check if not a tif ext?

    driver = gdal.GetDriverByName(drivercode)
    driver.Register()
    image = driver.Create(os.path.join(imagepath, imagename), 100, 100, 1, GDT_Int16)
    imageband = image.GetRasterBand(1)
    outarray = numpy.zeros(shape=(100, 50))
    outarray[outarray == 0] = 1000
    imageband.WriteArray(outarray)
    imageband = ""
    image = ""

    return 0


def get_hdf_subdatasets(hdfpath):
    #TODO docstrings

    hdf = gdal.Open(hdfpath, GA_ReadOnly)

    if hdf is None:
        raise Exception("Could not open " + hdfpath)

    subdatasets = []
    hdfsds = hdf.GetSubDatasets()

    for data in hdfsds:
        subdatasets.append((data[0], data[0].split(" ")[-1]))

    hdf = ""

    return subdatasets


def build_multiband_image(rootDIR, outName, newfoldername, find, drivercode, ndvalue, outputdir=None):
    """
    ##Set Args##
    rootdirectory = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2012/"
    outputfilename = "test"
    newfoldername = "kansas"
    VItofind = "EVI"
    drivercode = "ENVI"
    nodatavalue = -3000
    #projection = "PROJCS[\"Sinusoidal\",GEOGCS[\"GCS_Undefined\",DATUM[\"D_Undefined\",
                   SPHEROID[\"User_Defined_Spheroid\",6371007.181,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",
                   0.017453292519943295]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"False_Easting\",0.0],
                   PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",0.0],UNIT[\"Meter\",1.0]]"

    sys.exit(build_multiband_image(rootdirectory, outputfilename, newfoldername, VItofind, drivercode, nodatavalue)
    """

    #TODO docstrings

    if outputdir is None:
        outputdir = rootDIR

    outdir = create_output_dir(outputdir, newfoldername)
    print "\nOutputting files to : {0}".format(outdir)

    print "\nFinding HDF files in directory/subfolders: {0}".format(rootDIR)
    hdfs = find_files(rootDIR, ".hdf")
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

    #print "\nGetting output parameters..."
    #rows, cols, datatype, geotransform, projection = open_image(toprocess[0])
    #print "\tParameters: rows: {0}, cols: {1}, datatype: {2}, projection: {3}.".format(rows, cols, datatype, projection)

    outfile = os.path.join(outdir, outName)
    print "\nOutput file is: {0}".format(outfile)

    ## Create output file from first file to process ##
    template = openImage(toprocess[0])
    templateproperties = gdalProperties(template)
    outds = copySchemaToNewImage(templateproperties, outfile, numberofbands=bands, drivername=drivercode)
    template = ""
    del template
    print "\tCreated output file."

    print"\nAdding bands to output file..."
    for i in range(0, bands):
        print "\tProcessing band {0} of {1}...".format(i + 1, bands)
        print toprocess[i]
        image = openImage(toprocess[i])
        band = image.GetRasterBand(1)

        outband = outds.GetRasterBand(i + 1)

        print "\t\tReading band data to array..."
        data = band.ReadAsArray(0, 0, templateproperties.cols, templateproperties.rows)

        print "\t\tWriting band data to output band..."
        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(ndvalue)
        outband.FlushCache()

        outband = ""
        del data, outband
        band = ""
        image = ""

    print "\tFinished adding bands to output file."

    outds = ""
    del outds

    print "\nProcess completed."