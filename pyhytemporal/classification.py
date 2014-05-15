from osgeo import gdal
from osgeo.gdalconst import *

gdal.UseExceptions()


class GDAL_Object(object):
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
            self.geotransform = self.gdal.GetGeoTransform()
            self.projection = self.gdal.GetProjection()

            band = ""
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
            self.gdal = gdal.Open(infilepath, GA_ReadOnly)

            if self.gdal is None:
                raise Exception("Error encountered opening file.")
            else:
                self.hasParameters = True
                self.updateattributes()

        else:
            raise Exception("A file has already been opened with is object instance.")

        return

    def createNewImage(self, outfilepath, cols, rows, bands, datatype, drivername="GTiff", geotransform=None, projection=None):
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
                    self.gdal.SetGeotransform(geotransform)

                if projection:
                    self.gdal.SetProjection(projection)

                self.hasParameters = True
                self.updateattributes()

        return

    def copySchemaToNewImage(self, outfilepath, numberofbands=None, drivername="GTiff"):
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

        newimage = None

        if not self.hasParameters:
            raise Exception("Object instance has been created but no image has been opened.")
        else:

            if not numberofbands:
                numberofbands = self.bands

            newimage = GDAL_Object
            newimage.createNewImage(
                outfilepath,
                self.cols,
                self.rows,
                numberofbands,
                self.datatype,
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

class Temporal_Signature(object):
    def __init__(self, reffilepath):
        pass
