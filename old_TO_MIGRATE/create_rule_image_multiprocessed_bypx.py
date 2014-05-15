__author__ = 'phoetrymaster'

from osgeo import gdal
from osgeo.gdalconst import *
import os
import numpy
from scipy import interpolate
from scipy import optimize
import multiprocessing
from datetime import datetime as dt
import sys

gdal.UseExceptions()


########## METHODS ##########


def read_reference_file(filepath):
    error = 0

    print "\nReading reference file {0}".format(filepath)
    with open(filepath, "r") as f:
        values = {}
        try:
            for line in f:
                if not line.startswith("/"):
                    if len(line) >= 2:
                        print line
                        doy, vivalue = line.split(" ")
                        values[int(doy)] = float(vivalue)
                    else:
                        #line is not properly formatted
                        pass

        except ValueError:
            print "ERROR: File is not formatted correctly."
            error = 1
        except Exception as e:
            print e
            print "ERROR: Unknown problem reading reference file {0}".format(filepath)
            error = 2

    return error, values


##################################################
#
# Name, parameters, creator, date, major changes, algorithm
#
##################################################

def build_reference_dictionary(reffilelist):

    if (__debug__):
        print "ENTER build_reference_dictionary", reffilelist

    print "\nBuilding reference dictionary..."
    refs = {}
    for f in reffilelist:
        refname, ext = os.path.splitext(os.path.basename(f))
        errcode, values = read_reference_file(file)
        if not errcode:
            refs[refname] = values

    if (__debug__):
        print "EXIT  build_reference_dictionary", refs

    return refs




def create_output_raster(outFile, cols, rows, bands, datatype, drivername="GTiff"):
    driver = gdal.GetDriverByName(drivername)
    driver.Register()

    #Create Raster
    outds = driver.Create(outFile, cols, rows, bands, datatype)

    return outds


def create_output_dir(root, directoryname):
    dirpath = os.path.join(root, directoryname)

    if os.path.isdir(dirpath):
        count = 1
        dirpath_ = dirpath + "_"
        while True:
            dirpath = dirpath_ + str(count)
            count += 1
            if not os.path.isdir(dirpath):
                break

    os.makedirs(dirpath)
    return dirpath


def get_sort_dates_values(vals, threshold=-3000):
    """Gets the DOY dates (the keys) in a list from dictionary values and sorts those, placing them in chronological order
    (list x0). Then the function iterates over these values and gets the corresponding values, thresholding values if
    they are lower than an optional threshold value (-3000 default == nodata in MODIS imagery), then appending them to
    the list y. x and y are then returned."""

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
    return (1.0 / len(valsf) * numpy.sum(
        (valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve))) ** 2 for i in x0)) ** (
                    1.0 / 2.0)

def geometric(x, x0, valsf, interpolatedreferencecurve):
    return (numpy.product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve))) ** 2 for i in x0]) ** (
                    1.0 / 2)) ** (1 / float(len(valsf)))


def find_fit(valsf, interpolatedreferencecurve, bestguess, threshold, mthd, bnds="", mean=arithmetic):

    x0, y0 = get_sort_dates_values(valsf, threshold=threshold)

    #fun = lambda x: (product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve))) ** 2 for i in x0]) ** (1.0 / 2)) ** (1 / float(len(valsf)))

    res = optimize.minimize(arithmetic, (1, 1, bestguess), args=(x0, valsf, interpolatedreferencecurve), method=mthd, bounds=bnds)

    return res.fun, res.x, res.message


def interpolate_ref(refdict):

    x1, y1 = get_sort_dates_values(refdict)
    interpolated_curve = interpolate.splrep(x1, y1)

    return interpolated_curve


def process_pixel(bands, bestguess, col, cropname, doyinterval, fitmthd, img, interpolatedCurve, outarray, row,
                  startDOY, thresh, ndvalue, toprint, meantype=arithmetic):
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
        res, transforms, message = find_fit(valsf, interpolatedCurve, bestguess, thresh, fitmthd, bnds=bnds, mean=meantype)
        if toprint:
            print "\tPixel r{0}, c{1}: {2}: {3}, {4}, {5}".format(row, col, cropname, res, transforms, message)
    else:
        if toprint:
            "\tPixel r{0}, c{1}: {2}: NO DATA.".format(row, col, cropname)
        res = ndvalue

    outarray[row, col] = res

    return outarray


def process_reference(outputdir, cropname, refvalues, img, cols, rows, bands, geotransform, projection, drivercode, startDOY, doyinterval, bestguess, thresh, fitmthd, toprint, subset, ndvalue, meantype):

    try:
        #Create output rasters for each crop type to hold residual values from fit and arrays
        print "Creating {0} output raster...".format(cropname)
        outfileName = os.path.join(outputdir, cropname) + ".tif"
        outfile = create_output_raster(outfileName, cols, rows, 1, GDT_Float32, drivername=drivercode)
        outfile.SetGeoTransform(geotransform)
        outfile.SetProjection(projection)
        outdataset = outfile.GetRasterBand(1)
        outarray = numpy.zeros(shape=(rows, cols))
        outarray[outarray == 0] = ndvalue
        print "Created {0} raster.".format(cropname)


        #Interpolate reference values and return dict with key as type and curve as value
        interpolatedCurve = interpolate_ref(refvalues)


        #Iterate through each pixel and calculate the fit for each ref curve; write RMSE to array
        if subset:
            for col, row in subset:
                outarray = process_pixel(bands, bestguess, col, cropname, doyinterval, fitmthd, img, interpolatedCurve, outarray,
                                row, startDOY, thresh, ndvalue, toprint, meantype=meantype)
        else:
            for row in range(0, rows):
                for col in range(0, cols):
                    outarray = process_pixel(bands, bestguess, col, cropname, doyinterval, fitmthd, img, interpolatedCurve, outarray,
                                row, startDOY, thresh, ndvalue, toprint, meantype=meantype)



        #Write output array values to file
        print "Writing {0} output file...".format(cropname)
        outdataset.WriteArray(outarray, 0, 0)

        print "\nProcessing {0} finished.".format(cropname)

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


def phenological_classificaion(imagetoprocess, outputdirectory, createfoldername, referencevalues, gdaldrivercode, startDOY, doyinterval, threshold, bestguess, fitmethod, toprint=True, subset=[], meantype=arithmetic):
    start = dt.now()
    print start
    try:
        print "\nProcessing {0}...".format(imagetoprocess)
        outdir = create_output_dir(outputdirectory, createfoldername)
        print "Outputting files to {0} in {1}\n".format(createfoldername, outputdirectory)

        gdal.AllRegister()

        #Open multi-date image to analyze
        img = gdal.Open(imagetoprocess, GA_ReadOnly)

        if img is None:
            raise Exception("Could not open: {0}".format(imagetoprocess))

        #Get image properties
        cols  = img.RasterYSize
        rows  = img.RasterXSize
        bands = img.RasterCount
        geotransform = img.GetGeoTransform()
        projection = img.GetProjection()
        band = img.GetRasterBand(1)
        ndvalue = band.GetNoDataValue()

        print "Input image dimensions are {0} columns by {1} rows and contains {2} bands.".format(cols, rows, bands)

        processes = []
        for key, val in referencevalues.items():
            p = multiprocessing.Process(target=process_reference, args=(outdir, key, val, img, cols, rows, bands, geotransform, projection, gdaldrivercode, startDOY, doyinterval, bestguess, threshold, fitmethod, toprint, subset, ndvalue, meantype))
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

########## PROCEDURE ##########

if __name__ == '__main__':

    imagepath = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/ARC_Testing/ClipTesting/ENVI_1/test_clip_envi_3.dat"
    outdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/OutImages/"

    newfoldername = "Testing"

    drivercode = 'ENVI'
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

    sys.exit(phenological_classificaion(imagepath, outdir, newfoldername, refs, drivercode, startDOY, interval, threshold, bestguess, fitmthd, meantype=mean))