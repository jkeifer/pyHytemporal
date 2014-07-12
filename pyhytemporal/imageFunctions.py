import numpy
from osgeo import gdal
from osgeo.gdalconst import *
from core import gdalObject, gdalProperties
from utils import change_geotransform


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

    return array


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