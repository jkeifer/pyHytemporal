from __future__ import print_function
from datetime import datetime as dt
import multiprocessing
import os
import sys
import numpy
from scipy import interpolate, optimize
from core import gdalProperties
from imageFunctions import copySchemaToNewImage, openImage, read_image_into_array
from utils import band_number_to_doy, log
from vectorFunctions import get_px_coords_from_shapefile
from constants import *


#Logging level
LOGGINGLEVEL = INFO


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
    return (1.0 / len(valsf) * numpy.sum(
        (valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
        ** 2 for i in x0)) ** (1.0 / 2.0)


def geometric(x, x0, valsf, interpolatedreferencecurve):
    #TODO docstrings
    return (numpy.product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
                           ** 2 for i in x0]) ** (1.0 / 2)) ** (1 / float(len(valsf)))


def find_fit(valsf, interpolatedreferencecurve, bestguess, fitmethod=None, bounds=None, mean=None, threshold=None):
    #TODO docstrings

    if mean is None:
        mean = arithmetic
    else:
        mean = eval(mean)

    if fitmethod is None:
        fitmethod = 'SLSQP'

    x0, y0 = get_sort_dates_values(valsf, threshold=threshold)

    #fun = lambda x: (product([(valsf[i] - (x[0] * interpolate.splev((x[1] * (i + x[2])), interpolatedreferencecurve)))
    #                 ** 2 for i in x0]) ** (1.0 / 2)) ** (1 / float(len(valsf)))

    res = optimize.minimize(mean, (1, 1, bestguess), args=(x0, valsf, interpolatedreferencecurve),
                            method=fitmethod, bounds=bounds)

    return res.fun, res.x, res.message


def process_pixel(bestguess, col, cropname, doyinterval, fitmthd, array, interpolatedCurve, outarray, row,
                  startDOY, ndvalue, bounds, meantype=None, thresh=None, logging=None):
    #TODO docstrings

    logf = print

    if logging:
        try:
            logf = logging.log
        except:
            pass


    valsf = {}
    hasdata = True

    pixel = array[row, col]

    for i in range(pixel.size):
        measured = pixel[i]

        if measured == ndvalue:
            hasdata = False

        doy = band_number_to_doy(i+1, startDOY, doyinterval)
        valsf[doy] = measured

    if hasdata:
        res, transforms, message = find_fit(valsf, interpolatedCurve, bestguess, fitmthd, bounds=bounds, mean=meantype,
                                            threshold=thresh)

        logf("\tPixel r{0}, c{1}: {2}: {3}, {4}, {5}".format(row, col, cropname, res, transforms, message))

    else:
        logf("\tPixel r{0}, c{1}: {2}: NO DATA.".format(row, col, cropname))

        res = ndvalue

    outarray[row, col] = res

    return outarray


def process_reference(outputdir, signature, array, imageproperties, startDOY, doyinterval, bestguess, ndvalue, bounds,
                      meantype=None, subset=None, fitmthd=None, thresh=None):
    #TODO docstrings

    outfileName = os.path.join(outputdir, signature.name) + ".tif"
    logfilename = os.path.join(outputdir, signature.name) + ".txt"
    reflog = log(logfilename, verbosity=LOGGINGLEVEL)

    try:
        #Create output rasters for each crop type to hold residual values from fit and arrays
        reflog.log("Creating {0} output raster...".format(signature.name))
        outfileName = os.path.join(outputdir, signature.name) + ".tif"
        outfile = copySchemaToNewImage(imageproperties, outfileName, numberofbands=1)
        outdataset = outfile.GetRasterBand(1)
        outarray = numpy.zeros(shape=(imageproperties.rows, imageproperties.cols))
        outarray[outarray == 0] = ndvalue
        reflog.log("Created {0} raster.".format(signature.name))


        #Interpolate reference values and return dict with key as type and curve as value
        interpolatedCurve = interpolate.splrep(signature.daysofyear, signature.vivalues)

        #Iterate through each pixel and calculate the fit for each ref curve; write RMSE to array
        if subset:
            for row, col in subset:
                if row < 0 or col < 0 or row >= imageproperties.rows or col >= imageproperties.cols:
                    #print("Pixel not in processing extent; skipping.")
                    pass
                else:
                    outarray = process_pixel(bestguess, col, signature.name, doyinterval, fitmthd,
                                            array, interpolatedCurve, outarray,
                                            row, startDOY, ndvalue, bounds,
                                            meantype=meantype, thresh=thresh, logging=reflog)
        else:
            for row in range(0, imageproperties.rows):
                for col in range(0, imageproperties.cols):
                    outarray = process_pixel(bestguess, col, signature.name, doyinterval,
                                             fitmthd, array, interpolatedCurve, outarray,
                                             row, startDOY, ndvalue, bounds,
                                             meantype=meantype, thresh=thresh, logging=reflog)

        #Write output array values to file
        reflog.log("Writing {0} output file...".format(signature.name))
        outdataset.WriteArray(outarray, 0, 0)
        outdataset.SetNoDataValue(ndvalue)

        reflog.log("Processing {0} finished.".format(signature.name))

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        reflog.log(e)
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=reflog.file)
        try:
            print(row, col)
        except:
            pass

    finally:
        print("\nClosing files...")
        try:
            outdataset = None
        except:
            pass
        try:
            outfile = None
        except:
            pass


def fit_refs_to_image(imagetoprocess, outputdirectory, signaturecollection, startDOY,
                               doyinterval, bestguess, threshold=None, ndvalue=-3000, fitmethod=None, subset=None,
                               meantype=None, workers=4, timebounds=None, xbounds=None, ybounds=None):
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

    sys.exit(fit_refs_to_image(imagepath, outdir, newfoldername, refs, startDOY, interval, threshold, bestguess, fitmthd, meantype=mean))
    """
    #TODO docstrings

    start = dt.now()
    print(start)
    try:
        print("\nProcessing {0}...".format(imagetoprocess))
        print("Outputting files to {0}\n".format(outputdirectory))

        #Open multi-date image to analyze
        image = openImage(imagetoprocess)
        imageproperties = gdalProperties(image)

        print("Input image dimensions are {0} columns by {1} rows and contains {2} bands.".format(imageproperties.cols,
                                                                                                  imageproperties.rows,
                                                                                                  imageproperties.bands))

        array = read_image_into_array(image)  # Read all bands into a 3d array representing the image stack (x, y, time orientation)
        image = ""

        if timebounds:
            timebounds = (bestguess + timebounds[0], bestguess + timebounds[1])
        else:
            timebounds = (bestguess - 10, bestguess + 10)

        if not xbounds:
            xbounds = (0.6, 1.4)

        if not ybounds:
            ybounds = (0.6, 1.4)

        bounds = (xbounds, ybounds, timebounds)
        print(bounds)

        if subset:
            subset = get_px_coords_from_shapefile(imagetoprocess, subset)

        processes = []
        for signum, signature in enumerate(signaturecollection.signatures, start=1):
            p = multiprocessing.Process(target=process_reference,
                                        args=(outputdirectory, signature, array, imageproperties, startDOY, doyinterval,
                                              bestguess, ndvalue),
                                        kwargs={"subset": subset, "fitmthd": fitmethod, "meantype": meantype,
                                                "thresh": threshold, "bounds": bounds})

            #TODO: Problem with joining/starting processes--original thread closes before others are completed -- believe this is now fixed.

            p.start()
            processes.append(p)

            if len(processes) == workers:
                for p in processes:
                    p.join()
                    processes.remove(p)

        for p in processes:
                    p.join()

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print(e)
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)

    finally:
        print(dt.now() - start)