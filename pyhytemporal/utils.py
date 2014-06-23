import os
from math import floor
from osgeo import ogr
from core import *

#gdal.UseExceptions()


########## EXCEPTION CLASSES ##########


class ShapeDataError(Exception):
    """
    Error for wrong geometry type when loading shapefiles
    """
    pass


############## FUNCTIONS ##############


def find_files(searchdir, ext, recursive=True):
    foundfiles = []

    if recursive:
        for root, dirs, files in os.walk(searchdir):
            for f in files:
                if f.upper().endswith(ext.upper()):
                    foundfile = os.path.join(root, f)
                    foundfiles.append(foundfile)
    else:
        files = os.listdir(searchdir)
        for f in files:
            if f.upper().endswith(ext.upper()):
                foundfile = os.path.join(searchdir, f)
                foundfiles.append(foundfile)

    return foundfiles


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


def band_number_to_doy(bandnumber, startDOY, imageryinterval):
    calcdoy = (bandnumber - 1) * imageryinterval + startDOY

    if calcdoy > 0:
        calcdoy -= calcdoy % 365 % imageryinterval - 1

    return calcdoy


#############################################
#              IMAGE FUNCTIONS              #
#############################################


def openImage(infilepath):
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
    image = gdal.Open(infilepath, GA_ReadOnly)

    if image is None:
        raise Exception("Error encountered opening file.")

    return image


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
            raise Exception("Error encountered creating output file.")
        else:

            if geotransform:
                gdalImage.SetGeoTransform(geotransform)

            if projection:
                gdalImage.SetProjection(projection)

        return gdalImage


#############################################
#      GET PX COORDS FROM A SHAPEFILE       #
#############################################


def load_points(shapefile):
    """
    Returns a list of coordinate from points in an input shapefile.

    Required Argument(s):
        - shapefile: The path to a point-geometry shapefile

    Optional Argument(s):
        - None

    Returns:
        - points: A list of tuples with the x, y coords for each point in the input shapefile
    """
    #TODO Clean up commenting

    # Open shapeData
    shapeData = ogr.Open(validateShapePath(shapefile))
    # Validate shapeData
    validateShapeData(shapeData)
    # Get the first layer
    layer = shapeData.GetLayer()
    # Initialize
    points = []
    # For each point,
    for index in xrange(layer.GetFeatureCount()):
        # Get
        feature = layer.GetFeature(index)
        geometry = feature.GetGeometryRef()
        # Make sure that it is a point
        if geometry.GetGeometryType() != ogr.wkbPoint:
            raise ShapeDataError('This function only accepts point geometry.')
            # Get pointCoordinates
        pointCoordinates = geometry.GetX(), geometry.GetY()
        # Append
        points.append(pointCoordinates)
        # Cleanup
        feature.Destroy()
        # Cleanup
    shapeData.Destroy()
    # Return
    return points


def get_ref_from_shapefile(shapefile):
    """
    Gets a spatial reference from an input shapefile.

    Required Arguement(s):
        - shapefile: The path to a shapefile

    Optional Argument(s):
        - None

    Returns:
        - spatialref: The spatial reference of then input shapefile in proj4 format
    """

    # Open shapeData
    shapeData = ogr.Open(validateShapePath(shapefile))
    # Validate shapeData
    validateShapeData(shapeData)
    # Get the first layer
    layer = shapeData.GetLayer()
    # Get spatial reference as proj4
    spatialref = layer.GetSpatialRef()

    return spatialref


def getSpatialReferenceFromProj4(proj4):
    """Return GDAL spatial reference object from proj4 string"""
    #TODO: docstrings

    spatialReference = ogr.SpatialReference()
    spatialReference.ImportFromProj4(proj4)
    return spatialReference


# Validate

def validateShapePath(shapePath):
    """Validate shapefile extension"""
    #TODO: docstrings

    return os.path.splitext(str(shapePath))[0] + '.shp'


def validateShapeData(shapeData):
    """Make sure we can access the shapefile"""

    #TODO: docstrings

    # Make sure the shapefile exists
    if not shapeData:
        raise ShapeDataError('The shapefile is invalid')
        # Make sure there is exactly one layer
    if shapeData.GetLayerCount() != 1:
        raise ShapeDataError('The shapefile must have exactly one layer')


def check_spatial_refs(srs1, srs2):
    if srs1 == srs2:
        return False
    else:
        return True


def get_px_coords_from_points(raster, shapefile):
    """
    Takes geographic coordinates from a shapefile and finds the corresponding pixel coordinates on a raster.

    rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"
    #rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
    shp = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/SampleAreas/samplepoints2012_clip1_new.shp"

    print get_px_coords_from_points(rst, shp)
    """
    #TODO docstrings

    pointcoords = load_points(shapefile)
    ref1 = get_ref_from_shapefile(shapefile)

    image = gdalObject()
    image.open(raster)
    image.close()

    referror = check_spatial_refs(ref1, image.projection)

    if referror:
        print "WARNING: Spatial Reference of raster does not match points shapefile. Output may not be as expected. For best resutls ensure reference systems are identical."
        #TODO Change to use warnings module

    #Raster edge coords
    left = image.geotransform[0]
    top = image.geotransform[3]
    right = image.cols * image.geotransform[1] + image.geotransform[0]
    bottom = image.cols * image.geotransform[5] + image.geotransform[3]

    pxcoords = []
    for coords in pointcoords:
        x = int(floor(image.cols * (coords[0] - left) / (right - left)))
        y = int(floor(image.rows * (coords[1] - top) / (bottom - top)))
        pxcoords.append((x, y))

    return pxcoords


def change_geotransform(originaltransform, newxmin, newymin):
    """

    """
    #TODO Docstring

    newtransform = list(originaltransform)

    newtransform[0] += newtransform[1] * newxmin
    newtransform[3] += newtransform[5] * newymin

    return tuple(newtransform)


def clip_raster_to_extent(inraster, outraster, xmin, ymin, xextent, yextent):
    """

    """
    #TODO Docstring

    original = gdal.Open(inraster)
    original = gdalObject()
    original.open(inraster)

    newgeotransform = change_geotransform(original.geotransform, xmin, ymin)

    outds = original.copySchemaToNewImage(outraster, cols=xextent, rows=yextent, geotransform=newgeotransform)

    for i in range(1, original.bands + 1):
        band = original.gdal.GetRasterBand(i)
        outband = outds.gdal.GetRasterBand(i)

        nodatavalue = band.GetNoDataValue()
        data = band.ReadAsArray(xmin, ymin, xextent, yextent)

        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(nodatavalue)
        outband.FlushCache()

        del data, outband, band

    original.close()
    outds.close()

    return outraster