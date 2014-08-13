from math import floor
import os
from osgeo import ogr
from osgeo import osr
from core import gdalProperties, ShapeDataError

ogr.UseExceptions()


def getfieldindex(layer, fieldname):
    index = None
    layerdefn = layer.GetLayerDefn()

    for i in xrange(layerdefn.GetFieldCount()):
        if layerdefn.GetFieldDefn(i).GetName == fieldname:
            index = i
            break

    return index


def load_points(shapefile, outSpatialRef=None, fieldtoread=None):
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
    value = None  # to prevent reference before assignment

    if outSpatialRef:
        inSpatialRef = get_ref_from_shapefile(shapefile)
        coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # Open shapeData
    shapeData = ogr.Open(validateShapePath(shapefile))
    # Validate shapeData
    validateShapeData(shapeData)
    # Get the first layer
    layer = shapeData.GetLayer()

    # if fieldtoread:
    #     fieldindex = getfieldindex(fieldtoread)
    #     if fieldindex is None:
    #         raise ShapeDataError("fieldtoread was not found in shape data.")

    # Initialize
    points = []
    # For each point,
    for index in xrange(layer.GetFeatureCount()):
        # Get
        feature = layer.GetFeature(index)

        if fieldtoread:
            value = feature.GetField(fieldtoread)

        geometry = feature.GetGeometryRef()

        if coordTrans:
            geometry.Transform(coordTrans)

        # Make sure that it is a point
        if geometry.GetGeometryType() != ogr.wkbPoint:
            raise ShapeDataError('This function only accepts point geometry.')

        if value:
            pointCoordinates = geometry.GetX(), geometry.GetY(), value
        else:
            pointCoordinates = geometry.GetX(), geometry.GetY()

        points.append(pointCoordinates)
        # Cleanup
        feature.Destroy()
        # Cleanup
    shapeData.Destroy()
    # Return
    return points


def read_shapefile_to_points(shapefile, outSpatialRef=None):
    shapeData = ogr.Open(validateShapePath(shapefile))
    validateShapeData(shapeData)

    layer = shapeData.GetLayer()

    poly = layer.GetNextFeature()

    geom = poly.GetGeometryRef()

    if outSpatialRef:
        inSpatialRef = get_ref_from_shapefile(shapefile)
        coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        geom.Transform(coordTrans)

    extent = geom.GetEnvelope()
    points = geom.GetGeometryRef(0)

    return extent, points.Clone()


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

    return osr.SpatialReference().ImportFromProj4(proj4)


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


def get_px_coords_from_geographic_coords(gdalPropertiesObject, pointcoords):
    """

    """
    #TODO docstrings

    image = gdalPropertiesObject

    # get raster edge coords
    left = image.geotransform[0]
    top = image.geotransform[3]
    right = image.cols * image.geotransform[1] + image.geotransform[0]
    bottom = image.rows * image.geotransform[5] + image.geotransform[3]

    # calc px coords for each set of point coords
    pxcoords = []
    for coords in pointcoords:
        col = int(floor(image.cols * (coords[0] - left) / (right - left)))
        row = int(floor(image.rows * (coords[1] - top) / (bottom - top)))
        pxcoords.append((row, col))

    return pxcoords


def get_geographic_coords_from_px_coords(gdalPropertiesObject, pxcoords):
    """
    uses (row, col) format
    """
    #TODO docstrings

    image = gdalPropertiesObject

    # get raster edge coords
    left = image.geotransform[0]
    top = image.geotransform[3]
    horz_px_size = image.geotransform[1]
    vert_px_size = image.geotransform[5]


    # calc geo coords for each set of px coords
    pointcoords = []
    for coord in pxcoords:
        x = left + coord[1] * horz_px_size
        y = top + coord[0] * vert_px_size
        pointcoords.append((x, y))

    return pointcoords


def get_px_coords_from_shapefile(raster, shapefile):
    """
    Takes geographic coordinates from a shapefile and finds the corresponding pixel coordinates on a raster.

    rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"
    #rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
    shp = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/SampleAreas/samplepoints2012_clip1_new.shp"

    print get_px_coords_from_shapefile(rst, shp)
    """
    #TODO docstrings

    from imageFunctions import openImage

    # open, close image file and get properties
    raster = openImage(raster)
    imageproperties = gdalProperties(raster)

    rasterwkt = raster.GetProjectionRef()

    oSRSop = osr.SpatialReference()
    oSRSop.ImportFromWkt(rasterwkt)
    raster = None

    shppoints = load_points(shapefile, oSRSop)

    # get pixel coords from point coords
    pxcoords = get_px_coords_from_geographic_coords(imageproperties, shppoints)

    return pxcoords