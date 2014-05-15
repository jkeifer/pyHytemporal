__author__ = 'phoetrymaster'

from osgeo import gdal
from osgeo.gdalconst import *
import os
import numpy
from numpy import sum
from scipy import interpolate
from scipy import optimize
import multiprocessing
from datetime import datetime as dt

gdal.UseExceptions()

imagepath = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/ARC_Testing/ClipTesting/ENVI_1/test_clip_envi_3.dat"
rootdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/OutImages/"

newfoldername = "Testing"

drivercode = 'ENVI'
ndvalue = -3000

startDOY = 1
thresh = 500
bestguess = 0
fitmthd = 'SLSQP'

soylocs = [(6002, 2143), (5944, 2102), (5746, 2183), (5998, 2171)]
cornlocs = [(5997, 2139), (5940, 2096), (6051, 2230), (5691, 1998)]
wheatlocs = [(5993, 2136), (5937, 2080), (5935, 2076), (5921, 2217)]
refstoget = {"soy": soylocs, "corn": cornlocs, "wheat": wheatlocs}

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


########## METHODS ##########

class threads (multiprocessing.Process):
    def __init__(self, threadID, valsf, valsh, bestguess, threshold, fitmethod, key, row, col, fit, outarrays):
        multiprocessing.Process.__init__(self)
        self.threadID = threadID
        self.valsf = valsf
        self.valsh = valsh
        self.bestguess = bestguess
        self.thresh = threshold
        self.fitmethod = fitmethod
        self.key = key
        self.row = row
        self.col = col
        self.fit = fit
        self.outarrays = outarrays
    def run(self):
        self.res, self.transforms, self.message = find_fit(self.valsf, self.valsh, self.bestguess, threshhold=self.thresh, mthd=self.fitmethod)
        self.outarrays[self.key][self.row, self.col] = self.res
        print "\t{0}: {1}, {2}, {3}".format(self.key, self.res, self.transforms, self.message)
        self.fit[self.res] = self.threadID


def create_output_raster(outFile, cols, rows, bands, datatype, drivername="GTiff"):
    driver = gdal.GetDriverByName(drivername)
    driver.Register()

    outds = driver.Create(outFile, cols, rows, bands, datatype)

    return outds


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


def get_sort_dates_values(vals, threshhold=-3000):
    """Gets the DOY dates (the keys) in a list from dictonary vals and sorts those, placing them in chronological order
    (list x0). Then the function iterates over these values and gets the corresponding values, thresholding values if
    they are lower than an optional threshhold value (-3000 default = NoData in MODIS imagery), then appending them to
    the list y. x and y are then returned."""

    x = vals.keys()
    x.sort()
    y = []

    for i in x:
        if vals[i] < threshhold:
            y.append(threshhold)
        else:
            y.append(vals[i])

    return x, y


def find_fit(valsf, tck, bestguess, threshhold, mthd="TNC"):

    x0, y0 = get_sort_dates_values(valsf, threshhold=threshhold)

    fun = lambda x: ((1 / 22.8125 * sum(
        (valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), tck))) ** 2 for i in x0)) ** (
                    1. / 2))

    bnds = ((0.6, 1.4), (0.6, 1.4), (-10, 10))

    res = optimize.minimize(fun, (1, 1, bestguess), method=mthd, bounds=bnds)

    return res.fun, res.x, res.message


def interpolate_ref(refdict):

    x1, y1 = get_sort_dates_values(refdict)
    interpolatedcurve = interpolate.splrep(x1, y1)

    return interpolatedcurve


def main_processing(outputdir, cropname, refvalues, img, cols, rows, bands, geotransform, projection, drivercode):

    try:
        #Create output rasters for each crop type to hold residual values from fit and arrays
        print "\nCreating output files..."
        outfileName = os.path.join(outputdir, cropname) + ".tif"
        outfile = create_output_raster(outfileName, cols, rows, 1, GDT_Float32, drivername=drivercode)
        outfile.SetGeoTransform(geotransform)
        outfile.SetProjection(projection)
        outdataset = outfile.GetRasterBand(1)
        outarray = numpy.zeros(shape=(rows, cols))
        print "\tCreated file."


        #Interpolate reference values and return dict with key as type and curve as value
        interpolatedCurve = interpolate_ref(refvalues)

        #Iterate through each pixel and calculate the fit for each ref curve; write residuals and best fit to rasters
        print "\nFinding pixel fits..."
        for row in range(0, rows):
            arrays = {}
            for i in range(0, bands):
                band = img.GetRasterBand(i+1)
                arrays[i+1] = band.ReadAsArray(0, row, cols, 1)
            for col in range(0, cols):
                valsf = {}
                #print "Pixel r:{0}, c:{1}:".format(row, col)
                for i in range(0, bands):
                    measured = arrays[i+1][0, col]
                    valsf[startDOY + i*16] = measured
                res, transforms, message = find_fit(valsf, interpolatedCurve, bestguess, threshhold=thresh, mthd=fitmthd)
                outarray[row, col] = res
                #print "\t{0}: {1}, {2}, {3}".format(key, res, transforms, message)

        #Write output array values to files
        print "Writing {0} output file...".format(cropname)
        outdataset.WriteArray(outarray, 0, 0)

        print "\nProcessing {0} finished.".format(cropname)

    except Exception as e:
        print e

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

########## PROCEDURE ##########
start = dt.now()
print start

try:
    outputdirectory = create_output_dir(rootdir, newfoldername)
    print "\nOutputting files to : {0}".format(outputdirectory)

    gdal.AllRegister()

    #Open multi-date image to analyze
    img = gdal.Open(imagepath, GA_ReadOnly)

    if img is None:
        raise Exception("Could not open: {0}".format(imagepath))

    #Get image properties
    cols = img.RasterYSize
    rows = img.RasterXSize
    bands = img.RasterCount
    geotransform = img.GetGeoTransform()
    projection = img.GetProjection()

    print "Input image dimensions are {0} columns by {1} rows and contains {2} bands.".format(cols, rows, bands)

    processes = []
    for key, val in refs.items():
        p = multiprocessing.Process(target=main_processing, args=(outputdirectory, key, val, img, cols, rows, bands, geotransform, projection, drivercode))
        p.start()
        processes.append(p)

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