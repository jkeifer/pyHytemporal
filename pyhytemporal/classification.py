import os
import sys
import multiprocessing
import numpy
from math import floor
from datetime import datetime as dt
from osgeo import gdal
from osgeo.gdalconst import *
from osgeo import ogr
from scipy import interpolate
from scipy import optimize
from pyhytemporal.utils import *

gdal.UseExceptions()


########## EXCEPTION CLASSES ##########


class ShapeDataError(Exception):
    """
    Error for wrong geometry type when loading shapefiles
    """
    pass


########## OBJECT CLASSES #############


class gdalObject(object):
    """

    """
    #TODO docstrings

    def __init__(self):
        self.hasParameters = False
        self.gdal = None
        self.cols = None
        self.rows = None
        self.bands = None
        self.datatype = None
        self.geotransform = None
        self.projection = None

    def updateAttributes(self):
        """
        Gets the attributes of the raster and updates the properties of the instance.

        Required Argument(s):
            - None

        Optional Argument(s):
            - None

        Returns:
            - None
        """
        if self.gdal:
            self.rows = self.gdal.RasterYSize
            self.cols = self.gdal.RasterXSize
            self.bands = self.gdal.RasterCount
            band = self.gdal.GetRasterBand(1)
            self.datatype = band.DataType
            del band
            self.geotransform = self.gdal.GetGeoTransform()
            self.projection = self.gdal.GetProjection()

        else:
            raise Exception("No image is currently open from which the attributes can be read.")

        return

    def open(self, infilepath):
        """
        Opens a raster file and updates the instance attributes.

        Required Argument(s):
            - infilepath: The path to the raster to be opened.

        Optional Argument(s):
            - None

        Returns:
            - None
        """

        if not self.hasParameters:
            gdal.AllRegister()
            self.gdal = gdal.Open(infilepath, GA_ReadOnly)

            if self.gdal is None:
                raise Exception("Error encountered opening file.")
            else:
                self.hasParameters = True
                self.updateAttributes()

        else:
            raise Exception("A file has already been opened with is object instance. Close with reset=True to reuse.")

        return

    def createNewImage(self, outfilepath, cols, rows, bands, datatype,
                       drivername="GTiff", geotransform=None, projection=None):

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

        if not self.hasParameters:
            driver = gdal.GetDriverByName(drivername)
            driver.Register()

            self.gdal = driver.Create(outfilepath, cols, rows, bands, datatype)

            if self.gdal is None:
                raise Exception("Error encountered creating output file.")
            else:

                if geotransform:
                    self.gdal.SetGeoTransform(geotransform)

                if projection:
                    self.gdal.SetProjection(projection)

                self.hasParameters = True
                self.updateAttributes()

        return

    def copySchemaToNewImage(self, outfilepath, numberofbands=None, drivername=None, datatype=None):
        """
        Creates a new image using the properties of an existing image which has been loaded into a GDAL_Object.

        Required Argument(s):
            - outfilepath: The output path for the new file.

        Optional Argument(s):
            - numberofbands: If this is not specified, it will use the number of bands in the input image.
            - drivername: This is the file type for the output specified by GDAL driver name. The default is GeoTIFF.

        Returns:
            - GDAL_Object for new image.
        """

        if drivername is None:
            drivername = "GTiff"

        if datatype is None:
            datatype = self.datatype

        newimage = None

        if not self.hasParameters:
            raise Exception("Object instance has been created but no image has been opened.")
        else:

            if not numberofbands:
                numberofbands = self.bands

            newimage = gdalObject()
            newimage.createNewImage(
                outfilepath,
                self.cols,
                self.rows,
                numberofbands,
                datatype,
                drivername=drivername,
                geotransform=self.geotransform,
                projection=self.projection)

        return newimage

    def close(self, reset=False):
        """
        Closes an open image, but retains the image parameters for further use.

        Required Argument(s):
            - None

        Optional Argument(s):
            - reset: A boolean value. If true, will reset all image parameters to None.

        Returns:
            - None
        """

        self.gdal = None

        if reset:
            self.hasParameters = False
            self.cols = None
            self.rows = None
            self.bands = None
            self.datatype = None
            self.geotransform = None
            self.projection = None

        return

    def info(self):
        """
        Return the properties of the object.

        Required Argument(s):
            - None

        Optional Argument(s):
            - None

        Returns:
            - properties: A dict of all the properties of the object
        """

        properties = {"GDAL object": self.gdal, "Rows": self.rows, "Cols": self.cols, "Number of bands": self.bands,
                      "GDAL Datatype": self.datatype, "Geotransform": self.geotransform, "Projection": self.projection}

        return properties


class signatureCollection(object):
    """
    An object representing a collection of temporal signature objects.

    Properties:
        - self.signatures: A list of signature objects.

    Methods:
        - self.add(): Adds a signature object to the collection.
        - self.remove(): Removes a signature object from the collection.
    """

    def __init__(self, viName=None):
        """
        Initialize the signature collection.

        Required Argument(s):
            - None

        Optional Argument(s):
            - viName: The name or type of vegetation index from which the signatures were developed.

        Returns:
            - Nothing
        """

        self.signatures = []
        self.viName = viName

    def add(self, reffilepath, signaturename=None):
        """
        Reads a signature file (.ref) to create a signature object. This object is then appended to the signatures
        property on the instance.

        Required Argument(s):
            - reffilepath: The file path to the reference file (.ref)

        Optional Argument(s):
            - signaturename: This is the name of the signature, i.e. corn. If this is not provided, the reference file's
                name without the extension will be used.

        Returns:
            - signaturename
            - newsig.values: A tuple of tuples representing the DOY and measurment for each of the samples in the
                signature file.
        """

        if not signaturename:
            signaturename = os.path.splitext(os.path.basename(reffilepath))[0]

        newsig = temporalSignature(reffilepath, signaturename)
        self.signatures.append(newsig)

        return signaturename, newsig.values

    def remove(self, signaturename=None, index=None):
        """
        Removes a signature from the list of signatures. If no arguments are given, nothing will happen. The user must
        specify EITHER the signature name to be removed, or the index of the signature in the signature list. If both
        are provided, the program will not remove any signatures unless the index provided matches the index found when
        searching the list. If multiple signatures in the signature list have the same name, the one with the lowest
        index (the first to occur) will be removed, regardless if that is the signature the user desires to remove.

        Keyword Arguments:
            - signaturename: The name of the signature to be removed.
            - index: the index of the signature in the signature list that is to be removed.

        Returns:
            - removed: A boolean indicating whether or not the list item was removed.
        """

        removed = False
        toremove = None

        if not signaturename is None:
            foundindex = [i for i, j in enumerate(self.signatures) if j.name == signaturename][0]
            print foundindex
            if index:
                print "1:", index
                if index == foundindex:
                    toremove = index
                    print "2:", toremove
                else:
                    toremove = None
                    print "3:", toremove
            else:
                toremove = foundindex
                print "4:", toremove
        elif not index is None:
            toremove = index
            print "5:", toremove

        if not toremove is None:
            print "6: removing"
            del self.signatures[toremove]
            removed = True
        else:
            pass

        return removed


class temporalSignature(object):
    """
    This is an object representing a temporal signature of some plant/material.

    Properties:
        - self.name: The name of the plant/material. String
        - self.daysofyear: A tuple of integers representing the days of the year (DOYs)that have measurements.
        - self.vivalues: A tuple of floats that are the VI measurements on the DOYs in self.daysofyear.
        - self.values: A tuple of tuples for each of the measurements. The format is ((DOY1, VI), (DOY2, VI)).

    Methods:
        - ___init__: Reads an input .ref file, extracts the DOY and VI values, and creates the tuples.
    """

    def __init__(self, reffilepath, signaturename):
        """
        Extracts the DOY and VI values to lists, turns the lists to tuples to make them immutable, then adds them as
        properties to the object along with the name of the plant/material, and a tupled zip of the DOY and VI tuples.

        Required Argument(s):
            - reffilepath: The path to the .ref file with the signature to be imported.
            - signaturename: The name of the plant/material.

        Optional Argument(s):
            - None

        Returns:
            - Nothing
        """

        self.name = signaturename

        daysofyear, vivalues = [], []
        with open(reffilepath, "r") as f:
                for line in f:
                    if not line.startswith("/"):
                        if len(line) >= 2:
                            doy, vivalue = line.split(" ")
                            daysofyear.append(int(doy))
                            vivalues.append(float(vivalue))
                        else:
                            #line is not properly formatted
                            #raise Exception("Reference file is not formatted properly.")
                            pass

        if sorted(daysofyear) != daysofyear:
            print(sorted(daysofyear))
            print(daysofyear)
            raise Exception("Error: dates in reference file are not listed sequentially.")

        self.daysofyear = tuple(daysofyear)
        self.vivalues = tuple(vivalues)
        self.values = tuple(zip(self.daysofyear, self.vivalues))


################ FUNCTIONS ##################


#############################################
#           BUILD MULTIBAND IMAGE           #
#############################################


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


def build_multiband_image(rootDIR, outName, newfoldername, find, drivercode, ndvalue):
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

    outdir = create_output_dir(rootDIR, newfoldername)
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

    outfile = os.path.join(outdir, outName) + ".tif"
    print "\nOutput file is: {0}".format(outfile)

    ## Create output file from first file to process ##
    template = gdalObject()
    template.open(toprocess[0])
    outds = template.copySchemaToNewImage(outfile, numberofbands=bands, drivername=drivercode)
    template.close()
    del template
    print "\tCreated output file."

    print"\nAdding bands to output file..."
    for i in range(0, bands):
        print "\tProcessing band {0} of {1}...".format(i + 1, bands)
        image = gdal.Open(toprocess[i])
        band = image.GetRasterBand(1)

        outband = outds.gdal.GetRasterBand(i + 1)

        print "\t\tReading band data to array..."
        data = band.ReadAsArray(0, 0, outds.cols, outds.rows)

        print "\t\tWriting band data to output band..."
        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(ndvalue)
        outband.FlushCache()

        del data, outband
        image = ""

    print "\tFinished adding bands to output file."

    outds.close()
    del outds

    print "\nProcess completed."


#############################################
#     CREATE RULE IMAGE MULTIPROCESSED      #
#############################################


def get_sort_dates_values(vals, threshold=None):
    """Gets the DOY dates (the keys) in a list from dictionary values and sorts those, placing them in chronological order
    (list x0). Then the function iterates over these values and gets the corresponding values, thresholding values if
    they are lower than an optional threshold value (-3000 default == nodata in MODIS imagery), then appending them to
    the list y. x and y are then returned."""
    #TODO docstrings

    if threshold is None:
        threshold = -3000

    x = vals.keys()
    x.sort()
    y = []

    for i in x:
        if vals[i] < threshold:
            y.append(threshold)
        else:
            y.append(vals[i])

    return x, y


def arithmetic(x, x0, valsf, interpolatedreferencecurve):
    #TODO docstrings
    return (1.0 / len(valsf) * numpy.sum((valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
            ** 2 for i in x0)) ** (1.0 / 2.0)


def geometric(x, x0, valsf, interpolatedreferencecurve):
    #TODO docstrings
    return (numpy.product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
                           ** 2 for i in x0]) ** (1.0 / 2)) ** (1 / float(len(valsf)))


def find_fit(valsf, interpolatedreferencecurve, bestguess, fitmethod=None, bounds=None, mean=None, threshold=None):
    #TODO docstrings

    if mean is None:
        mean = arithmetic

    if fitmethod is None:
        fitmethod = 'SLSQP'

    x0, y0 = get_sort_dates_values(valsf, threshold=threshold)

    #fun = lambda x: (product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
    #                 ** 2 for i in x0]) ** (1.0 / 2)) ** (1 / float(len(valsf)))

    res = optimize.minimize(mean, (1, 1, bestguess), args=(x0, valsf, interpolatedreferencecurve),
                            method=fitmethod, bounds=bounds)

    return res.fun, res.x, res.message


def process_pixel(bands, bestguess, col, cropname, doyinterval, fitmthd, img, interpolatedCurve, outarray, row,
                  startDOY, ndvalue, meantype=None, thresh=None):
    #TODO docstrings

    valsf = {}
    i2 = 0
    st = startDOY
    hasdata = True

    for i in range(0, bands):
        band = img.GetRasterBand(i + 1)
        measured = int(band.ReadAsArray(col, row, 1, 1))

        if measured == ndvalue:
            hasdata = False

        doy = st + i2 * doyinterval

        if  doy > 365:
            doy = 366
            i2  = 0
            st  = 366

        valsf[doy] = measured
        i2 += 1

    if hasdata:
        bnds = ((0.6, 1.4), (0.6, 1.4), (-10, 10))
        res, transforms, message = find_fit(valsf, interpolatedCurve, bestguess, fitmthd, bounds=bnds, mean=meantype, threshold=thresh)

        if __debug__:
            print "\tPixel r{0}, c{1}: {2}: {3}, {4}, {5}".format(row, col, cropname, res, transforms, message)

    else:

        if __debug__:
            "\tPixel r{0}, c{1}: {2}: NO DATA.".format(row, col, cropname)

        res = ndvalue

    outarray[row, col] = res

    return outarray


def process_reference(outputdir, signature, img, startDOY, doyinterval, bestguess, ndvalue,
                      meantype=None, drivercode=None, subset=None, fitmthd=None, thresh=None):
    #TODO docstrings

    try:
        #Create output rasters for each crop type to hold residual values from fit and arrays
        print "Creating {0} output raster...".format(signature.name)
        outfileName = os.path.join(outputdir, signature.name) + ".tif"
        outfile = img.copySchemaToNewImage(outfileName, numberofbands=1, drivername=drivercode)
        outdataset = outfile.gdal.GetRasterBand(1)
        outarray = numpy.zeros(shape=(outfile.rows, outfile.cols))
        outarray[outarray == 0] = ndvalue
        print "Created {0} raster.".format(signature.name)


        #Interpolate reference values and return dict with key as type and curve as value
        interpolatedCurve = interpolate.splrep(signature.daysofyear, signature.vivalues)

        #Iterate through each pixel and calculate the fit for each ref curve; write RMSE to array
        if subset:
            for col, row in subset:
                outarray = process_pixel(img.bands, bestguess, col, signature.name, doyinterval, fitmthd, img.gdal, interpolatedCurve, outarray,
                                row, startDOY, thresh, ndvalue, meantype=meantype)
        else:
            for row in range(0, img.rows):
                for col in range(0, img.cols):
                    outarray = process_pixel(img.bands, bestguess, col, signature.name, doyinterval, fitmthd, img.gdal, interpolatedCurve, outarray,
                                row, startDOY, thresh, ndvalue, meantype=meantype)



        #Write output array values to file
        print "Writing {0} output file...".format(signature.name)
        outdataset.WriteArray(outarray, 0, 0)

        print "\nProcessing {0} finished.".format(signature.name)

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print e
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)

    finally:
        print "\nClosing files..."
        try:
            outdataset = None
        except:
            pass
        try:
            outfile = None
        except:
            pass


def phenological_classificaion(imagetoprocess, outputdirectory, createfoldername, signaturecollection, startDOY,
                               doyinterval, bestguess, threshold=None, ndvalue=-3000, fitmethod=None, subset=None, gdaldrivercode=None,
                               meantype=None, workers=4):
    """
    imagepath = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/ARC_Testing/ClipTesting/ENVI_1/test_clip_envi_3.dat"
    outdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/OutImages/"

    newfoldername = "Testing"

    drivercode1 = 'ENVI'
    #ndvalue = -3000

    startDOY = 1
    interval = 16
    threshold = 500
    bestguess = 0
    fitmthd = 'SLSQP'
    mean = geometric  # Acceptable values are geometric (geometric mean) and arithmetic (arithmetic mean). It is an optional argument for the classifier.


    refs = {
        'soy': {1: 174.5, 97: 1252.25, 65: 1139.5, 209: 7659.0, 273: 4606.75, 337: 1371.75, 17: 1055.5, 33: 1098.0,
                49: 1355.25,
                129: 1784.75, 257: 6418.0, 321: 1644.5, 305: 1472.75, 193: 5119.75, 289: 1878.75, 177: 3439.5, 241: 7565.75,
                81: 1205.5, 225: 7729.75, 145: 1736.25, 161: 1708.25, 353: 1358.25, 113: 1340.0},
        'corn': {1: 392.25, 97: 1433.25, 65: 1258.5, 209: 6530.0, 273: 1982.5, 337: 1658.5, 17: 1179.25, 33: 1196.75,
                 49: 1441.25, 129: 1885.25, 257: 2490.25, 321: 1665.75, 305: 1439.0, 193: 6728.25, 289: 1634.5,
                 177: 6356.75,
                 241: 4827.25, 81: 1355.75, 225: 5547.5, 145: 2196.5, 161: 3143.25, 353: 1704.75, 113: 1716.5},
        'wheat': {1: 719.75, 97: 6594.75, 65: 1935.25, 209: 2013.5, 273: 1493.5, 337: 1498.25, 17: 1816.5, 33: 1815.0,
                  49: 1985.25, 129: 6758.0, 257: 1685.75, 321: 1582.5, 305: 1163.25, 193: 2186.25, 289: 1264.5, 177: 2222.5,
                  241: 2301.0, 81: 4070.5, 225: 1858.0, 145: 6228.5, 161: 3296.5, 353: 1372.5, 113: 7035.25}
    }

    sys.exit(phenological_classificaion(imagepath, outdir, newfoldername, refs, drivercode1, startDOY, interval, threshold, bestguess, fitmthd, meantype=mean))
    """
    #TODO docstrings

    start = dt.now()
    print start
    try:
        print "\nProcessing {0}...".format(imagetoprocess)
        outdir = create_output_dir(outputdirectory, createfoldername)
        print "Outputting files to {0} in {1}\n".format(createfoldername, outputdirectory)

       #Open multi-date image to analyze
        img = gdalObject
        img.open(imagetoprocess)

        print "Input image dimensions are {0} columns by {1} rows and contains {2} bands.".format(img.cols, img.rows, img.bands)

        if subset:
            subset = get_px_coords_from_points(img.gdal, subset)

        processes = []
        for signature in signaturecollection.signatures:
            p = multiprocessing.Process(target=process_reference,
                                        args=(outdir, signature, img, startDOY, doyinterval, bestguess, ndvalue),
                                        kwargs={drivercode: gdaldrivercode, subset: subset, fitmthd: fitmethod,
                                                meantype: meantype, thresh: threshold})
            p.start()
            processes.append(p)

            if len(processes) == workers:
                for p in processes:
                    p.join()

        print dt.now() - start

    except Exception as e:
        print e

    finally:
        print "\nClosing file..."
        try:
            img = None
        except:
            pass


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

    image = gdalObject
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


#############################################
#       Classify and Assess Accuracy        #
#############################################


def generate_thresholds(start, step, numberofsteps, lengthofelement):
    #TODO Docstring

    from itertools import product
    end = numberofsteps * step + start
    thresholdvals = range(start, end, step)
    thresholdlists = []
    for i in range(0, lengthofelement):
        thresholdlists.append(thresholdvals)
    for i in product(*thresholdlists):
        yield i


def classify_with_threshold(croparray, filelist, searchdir, searchstringsvals, thresh, nodata):

    #TODO test to ensure thresh length is equal to the number of image files?

    arrays = []
    #print filelist

    i = 0
    for f, cropval in filelist:
        img = gdal.Open(f, GA_ReadOnly)
        if img is None:
            raise Exception("Could not open: {0}".format(os.path.join(searchdir, f)))
        else:
            rows = img.RasterYSize
            cols = img.RasterXSize
            band = img.GetRasterBand(1)
            array = band.ReadAsArray(0, 0, cols, rows)
            array[array > thresh[i]] = 10000
            arrays.append((numpy.copy(array), cropval))
            band = ""
            img = ""
        i += 1

    count = 0
    finals = []

    for array, cropval in arrays:
        ltarrays = []
        nodataarrays = []
        for i in range(0, len(arrays)):
            if not i == count:
                lt = array.__lt__(arrays[i][0])
                ltarrays.append(numpy.copy(lt))
                ndarray = numpy.copy(array)
                ndarray[ndarray != nodata] = 0
                ndarray = numpy.rint(ndarray)
                ndarray = ndarray.astype(int)
                nodataarrays.append(ndarray)
        count += 1

        for i in range(0, len(ltarrays)):
            if not i:
                allpxbestfit = numpy.copy(ltarrays[i])
            else:
                allpxbestfit = allpxbestfit.__and__(ltarrays[i])
        finals.append(cropval * allpxbestfit)

    nodataarray = ""
    for ndarray in nodataarrays:
        if nodataarray == "":
            nodataarray = numpy.copy(ndarray)
        else:
            nodataarray = nodataarray.__and__(ndarray)

    classification = ""
    for final in finals:
        if classification == "":
            classification = numpy.copy(final)
        else:
            classification = classification.__or__(final)
    classification = classification.__or__(nodataarray)

    #Accuracy Assessment
    results = {}
    for string, val in searchstringsvals:
        searchdict = {}
        temparray = numpy.copy(classification)
        temparray[temparray != val] = 0
        correct = temparray.__eq__(croparray)
        incorrect = temparray.__ne__(croparray)
        temparray[temparray == val] = 1
        incorrectvals = incorrect.__mul__(croparray).__mul__(temparray)
        for string2, val2 in searchstringsvals:
            temparray2 = numpy.copy(incorrectvals)
            temparray2[temparray2 != val2] = 0
            searchdict[string2] = temparray2.sum() / val2

        searchdict[string] = correct.sum()
        searchdict["other"] = len([x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
        results[string] = searchdict.copy()

    searchdict = {}
    temparray = numpy.copy(classification)
    temparray[classification == 0] = 1
    temparray[classification != 0] = 0
    incorrectvals = temparray.__mul__(croparray)

    for string2, val2 in searchstringsvals:
        temparray2 = numpy.copy(incorrectvals)
        temparray2[temparray2 != val2] = 0
        searchdict[string2] = temparray2.sum() / val2

    searchdict["other"] = len([x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
    results["other"] = searchdict.copy()

    numpx = 0
    for key, val in results.items():
        for k, v in val.items():
            numpx += v

    printstring = ""
    correct = 0
    croporder = []
    h = 0
    for crop, values in results.items():
        total = 0
        vals = []
        for crop2, pxcount in values.items():
            total = total + pxcount
            if not crop2 in croporder:
                croporder.append(crop2)
            if crop2 == crop:
                correct += pxcount
            vals.append(pxcount)
        printstring = printstring + "{0}\t\t{1}\t{2}\n".format(crop, vals, total)
        h += total
    accuracy = correct / (h * 1.0)
    outstring = ("{0}\n\t\t{1}\trow total\n{2}\n{3}\n\n\n".format(thresh, croporder, printstring, accuracy))

    return accuracy, classification, cols, rows, outstring


def main(searchdir, cropimgpath, searchstringsvals, nodata, outdir=None, outfilename=None):

    if outdir is None:
        outdir = searchdir
    else:
        pass

    if outfilename is None:
        today = dt.now()
        outfilename = today.strftime("%Y-%m-%d_%H%M%S_") + os.path.splitext(os.path.basename(cropimgpath))[0]


    outFile = os.path.join(outdir, outfilename + ".tif")

    accuracyreport = os.path.join(outdir, outfilename + ".txt")

    try:
        gdal.AllRegister()
        #np.set_printoptions(threshold=np.nan)

        #Crop image is constant for all iterations
        cropimg = gdalObject
        cropimg.open(cropimgpath)
        band = cropimg.gdal.GetRasterBand(1)
        croparray = band.ReadAsArray(0, 0, cropimg.cols, cropimg.rows)
        band = None
        cropimg.close()
        print("Opened crop img")

        filelist = []
        files = os.listdir(searchdir)
        for f in files:
            for string, val in searchstringsvals:
                if f.endswith(".tif"):
                    if string in f:
                        filelist.append((os.path.join(searchdir, f), val))

        thresholds = generate_thresholds(600, 100, 10, len(filelist))

        writestring = ""
        bestacc = 0
        bestthresh = ""
        for thresh in thresholds:
            start = dt.now()
            accuracy, classification, cols, rows, outstring = classify_with_threshold(croparray,
                                                                                           filelist,
                                                                                           searchdir, searchstringsvals,
                                                                                           thresh, nodata)
            writestring = writestring + outstring

            if accuracy > bestacc:
                bestacc = accuracy
                bestthresh = thresh

            elapsed = dt.now() - start
            print thresh, elapsed, accuracy, bestacc

    except Exception as e:
        print e

    finally:

        with open(accuracyreport, 'w') as text:
            text.write("Classification using curves from {0}".format(searchdir))
            text.write("{0}\nBest:\n{1} {2}".format(writestring, bestthresh, bestacc))

        print "\n", bestthresh, bestacc

        accuracy, classification, cols, rows, outstring = classify_with_threshold(croparray, filelist, searchdir,
                                                                                  searchstringsvals, bestthresh,
                                                                                  nodata)

        driver = gdal.GetDriverByName("ENVI")
        driver.Register()

        outds = driver.Create(outFile, cols, rows, 1, GDT_Int16)
        outds.SetGeoTransform(cropimg.geotransform)
        outds.SetProjection(cropimg.projection)
        outband = outds.GetRasterBand(1)
        outband.WriteArray(classification, 0, 0)
        outband.SetNoDataValue(nodata)
        outband.FlushCache()

        outband = ""
        outds = ""

        print "outputted"

        return 0


#############################################
#         GET TEMPORAL SIGNATURES           #
#############################################


def get_crop_pixel_values(imagepath, locations):

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
        for i in range(0, bands):
            band = img.GetRasterBand(i + 1)
            v = int(band.ReadAsArray(int(floor(location[0])), int(floor(location[1])), 1, 1))
            print "\tBand {0}: {1}".format(i + 1, v)
            values.append(v)
            band = None
        refs.append(values)
    img = None

    return refs


def write_refs_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing pixel curves to output file:"
    print referencevalues

    #output individual pixel curves to file
    with open(os.path.join(outdir, cropname + postfix + "_points.ref"), "w") as f:
        if comment:
            f.write("//"+comment+"\n\n")
        point = 1
        for points in referencevalues:
            f.write("\nPoint {0}:\n".format(point))
            point += 1
            imgnumber = 0
            st = startdoy
            for val in points:
                doy = st + imgnumber * doyinterval
                if doy > 365:
                    doy = 366
                    st = 366
                    imgnumber = 0
                f.write("{0} {1}\n".format(doy, val))
                imgnumber += 1


def write_mean_ref_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing mean reference curve to output file:"
    print referencevalues

    #output mean pixel values to file
    with open(os.path.join(outdir, cropname + postfix + "_mean.ref"), "w") as f:
        if comment:
            f.write("//"+comment+"\n\n")
        meanvals = get_mean_values(referencevalues)
        imgnumber = 0
        st = startdoy
        for val in meanvals:
            doy = int(startdoy) + imgnumber * int(doyinterval)
            if doy > 365:
                doy = 366
                st = 366
                imgnumber = 0
            f.write("{0} {1}\n".format(doy, val))
            imgnumber += 1


def get_mean_values(referencevalues):
    array = numpy.array(referencevalues)
    mean = numpy.mean(array, axis=0)
    return list(mean)


def get_reference_curves(image, refstoget, startdoy, imageinterval, outdir="", filepostfix=""):
    if not outdir:
        outdir = os.path.dirname(image)

    for key, val in refstoget.items():
        if val:
            cropname = key
            locs = val
            refvals = get_crop_pixel_values(image, locs)
            comment = "Generated from {0} by get_crop_pixel_values version 0.1.".format(image)
            write_refs_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)
            write_mean_ref_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)

########## PROCEDURE ##########

if __name__ == '__main__':
    sys.exit()