from osgeo import ogr
import os, sys
from build_multiband_image import open_image
from math import floor


def load_points(shapefile):
    """Given a shapePath, return a list of points in GIS coordinates"""

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
            raise ShapeDataError('This module can only load points; use geometry_store.py')
        # Get pointCoordinates
        pointCoordinates = geometry.GetX(), geometry.GetY()
        # Append
        points.append(pointCoordinates)
        # Cleanup
        feature.Destroy()
    # Get spatial reference as proj4
    spatialref = layer.GetSpatialRef()
    # Cleanup
    shapeData.Destroy()
    # Return
    return points, spatialref


def getSpatialReferenceFromProj4(proj4):
    """Return GDAL spatial reference object from proj4 string"""

    spatialReference = osr.SpatialReference()
    spatialReference.ImportFromProj4(proj4)
    return spatialReference


# Validate

def validateShapePath(shapePath):
    """Validate shapefile extension"""

    return os.path.splitext(str(shapePath))[0] + '.shp'


def validateShapeData(shapeData):
    """Make sure we can access the shapefile"""

    # Make sure the shapefile exists
    if not shapeData:
        raise ShapeDataError('The shapefile is invalid')
    # Make sure there is exactly one layer
    if shapeData.GetLayerCount() != 1:
        raise ShapeDataError('The shapefile must have exactly one layer')


def check_spatial_refs(srs1, srs2):
    if srs1 == srs2:
        return 0
    else:
        return 1


# Error

class ShapeDataError(Exception):
    pass


def get_px_coords_from_points(raster, shapefile):
    """Takes geographic coordinates from a shapefile and finds the corresponding pixel coordinates on a raster."""

    pointcoords, ref1 = load_points(shapefile)
    rows, cols, bandtype, geotransform, ref2 = open_image(raster)
    referror = check_spatial_refs(ref1, ref2)

    if referror:
        print "WARNING: Spatial Reference of raster does not match points shapefile. Output may not be as expected. For best resutls ensure reference systems are identical."

    #Raster edge coords
    left = geotransform[0]
    top = geotransform[3]
    right = cols * geotransform[1] + geotransform[0]
    bottom = cols * geotransform[5] + geotransform[3]

    pxcoords = []
    for coords in pointcoords:
        x = int(floor(cols * (coords[0] - left) / (right - left)))
        y = int(floor(rows * (coords[1] - top) / (bottom - top)))
        pxcoords.append((x, y))

    return pxcoords


if __name__ == '__main__':
    rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"
    #rst = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
    shp = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/SampleAreas/samplepoints2012_clip1_new.shp"

    print get_px_coords_from_points(rst, shp)