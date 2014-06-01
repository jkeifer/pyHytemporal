import os
import sys
from osgeo import gdal
from osgeo.gdalconst import *

gdal.UseExceptions()

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

    def copySchemaToNewImage(self, outfilepath, numberofbands=None, drivername=None, datatype=None, cols=None, rows=None,
                             geotransform=None, projection=None):
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

        if cols is None:
            cols = self.cols

        if rows is None:
            rows = self.rows

        if geotransform is None:
            geotransform = self.geotransform

        if projection is None:
            projection = self.projection

        newimage = None

        if not self.hasParameters:
            raise Exception("Object instance has been created but no image has been opened.")
        else:

            if not numberofbands:
                numberofbands = self.bands

            newimage = gdalObject()
            newimage.createNewImage(
                outfilepath,
                cols,
                rows,
                numberofbands,
                datatype,
                drivername=drivername,
                geotransform=geotransform,
                projection=projection)

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


if __name__ == '__main__':
    sys.exit()
