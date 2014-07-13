import os
import sys
import multiprocessing
from math import floor
from datetime import datetime as dt

import numpy
from osgeo import gdal
from osgeo.gdalconst import *
from scipy import interpolate
from scipy import optimize

from imageFunctions import *
from vectorFunctions import get_px_coords_from_points
from utils import *
from core import *
from plotting import Plot

lock = multiprocessing.Lock()

#gdal.UseExceptions()

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


def build_multiband_image(rootDIR, outName, newfoldername, find, drivercode, ndvalue, outputdir=None):
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

    if outputdir is None:
        outputdir = rootDIR

    outdir = create_output_dir(outputdir, newfoldername)
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

    outfile = os.path.join(outdir, outName)
    print "\nOutput file is: {0}".format(outfile)

    ## Create output file from first file to process ##
    template = openImage(toprocess[0])
    templateproperties = gdalProperties(template)
    outds = copySchemaToNewImage(templateproperties, outfile, numberofbands=bands, drivername=drivercode)
    template = ""
    del template
    print "\tCreated output file."

    print"\nAdding bands to output file..."
    for i in range(0, bands):
        print "\tProcessing band {0} of {1}...".format(i + 1, bands)
        print toprocess[i]
        image = openImage(toprocess[i])
        band = image.GetRasterBand(1)

        outband = outds.GetRasterBand(i + 1)

        print "\t\tReading band data to array..."
        data = band.ReadAsArray(0, 0, templateproperties.cols, templateproperties.rows)

        print "\t\tWriting band data to output band..."
        outband.WriteArray(data, 0, 0)
        outband.SetNoDataValue(ndvalue)
        outband.FlushCache()

        outband = ""
        del data, outband
        band = ""
        image = ""

    print "\tFinished adding bands to output file."

    outds = ""
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
                  startDOY, ndvalue, meantype=None, thresh=None):
    #TODO docstrings

    valsf = {}
    hasdata = True

    pixel = array[row, col]
    #print pixel
    #print pixel.size

    for i in range(pixel.size):
        measured = pixel[i]
        #print(measured, col, row, i + 1)

        if measured == ndvalue:
            hasdata = False

        doy = band_number_to_doy(i+1, startDOY, doyinterval)
        valsf[doy] = measured

    if hasdata:
        bnds = ((0.6, 1.4), (0.6, 1.4), (-10, 10))
        res, transforms, message = find_fit(valsf, interpolatedCurve, bestguess, fitmthd, bounds=bnds, mean=meantype,
                                            threshold=thresh)

        if __debug__:
            print "\tPixel r{0}, c{1}: {2}: {3}, {4}, {5}".format(row, col, cropname, res, transforms, message)

    else:

        if __debug__:
            "\tPixel r{0}, c{1}: {2}: NO DATA.".format(row, col, cropname)

        res = ndvalue

    outarray[row, col] = res

    return outarray


def process_reference(outputdir, signature, array, imageproperties, startDOY, doyinterval, bestguess, ndvalue,
                      meantype=None, subset=None, fitmthd=None, thresh=None):
    #TODO docstrings

    try:
        #Create output rasters for each crop type to hold residual values from fit and arrays
        print "Creating {0} output raster...".format(signature.name)
        outfileName = os.path.join(outputdir, signature.name) + ".tif"
        outfile = copySchemaToNewImage(imageproperties, outfileName, numberofbands=1)
        outdataset = outfile.GetRasterBand(1)
        outarray = numpy.zeros(shape=(imageproperties.rows, imageproperties.cols))
        outarray[outarray == 0] = ndvalue
        print "Created {0} raster.".format(signature.name)


        #Interpolate reference values and return dict with key as type and curve as value
        interpolatedCurve = interpolate.splrep(signature.daysofyear, signature.vivalues)

        #Iterate through each pixel and calculate the fit for each ref curve; write RMSE to array
        if subset:
            for row, col in subset:
                outarray = process_pixel(bestguess, col, signature.name, doyinterval, fitmthd,
                                         array, interpolatedCurve, outarray,
                                         row, startDOY, ndvalue, meantype=meantype, thresh=thresh)
        else:
            for row in range(0, imageproperties.rows):
                for col in range(0, imageproperties.cols):
                    outarray = process_pixel(bestguess, col, signature.name, doyinterval,
                                             fitmthd, array, interpolatedCurve, outarray,
                                             row, startDOY, ndvalue, meantype=meantype, thresh=thresh)

        #Write output array values to file
        print "Writing {0} output file...".format(signature.name)
        outdataset.WriteArray(outarray, 0, 0)
        outdataset.SetNoDataValue(ndvalue)

        print "\nProcessing {0} finished.".format(signature.name)

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print e
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)
        try:
            print(row, col)
        except:
            pass

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


def phenological_classificaion(imagetoprocess, outputdirectory, signaturecollection, startDOY,
                               doyinterval, bestguess, threshold=None, ndvalue=-3000, fitmethod=None, subset=None,
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

    sys.exit(phenological_classificaion(imagepath, outdir, newfoldername, refs, startDOY, interval, threshold, bestguess, fitmthd, meantype=mean))
    """
    #TODO docstrings

    start = dt.now()
    print start
    try:
        print "\nProcessing {0}...".format(imagetoprocess)
        print "Outputting files to {0}\n".format(outputdirectory)

        #Open multi-date image to analyze
        image = openImage(imagetoprocess)
        imageproperties = gdalProperties(image)

        print "Input image dimensions are {0} columns by {1} rows and contains {2} bands.".format(imageproperties.cols,
                                                                                                  imageproperties.rows,
                                                                                                  imageproperties.bands)

        array = read_image_into_array(image)  # Read all bands into a 3d array representing the image stack (time, x, y orientation)
        #array = array.transpose(1, 2, 0)  # Turn the array to allow a single pixel to be isolated in the stack (x, y, time orientation)

        if subset:
            subset = get_px_coords_from_points(imagetoprocess, subset)

        processes = []
        for signature in signaturecollection.signatures:
            p = multiprocessing.Process(target=process_reference,
                                        args=(outputdirectory, signature, array, imageproperties, startDOY, doyinterval,
                                              bestguess, ndvalue),
                                        kwargs={"subset": subset, "fitmthd": fitmethod, "meantype": meantype,
                                                "thresh": threshold})

            #TODO: Problem with joining/starting processes--original thread closes before others are completed

            p.start()
            processes.append(p)

            if len(processes) == workers:
                for p in processes:
                    p.join()
                    processes.remove(p)

        print dt.now() - start

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print e
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)

    finally:
        print dt.now() - start
        print "\nClosing file..."
        try:
            image = None
        except:
            pass


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

    i = 0  # Track the iterations to set the thresholds for each image
    for imagefile, cropval in filelist:
        img = openImage(imagefile)
        #TODO Convert to new properties object and refactor into smaller, testable methods...
        if img is None:
            raise Exception("Could not open: {0}".format(os.path.join(searchdir, imagefile)))
        else:
            band = img.GetRasterBand(1)
            array = band.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize)
            array[array > thresh[i]] = 10000  # Make all values in the array greater than the thresh equal 10000
            arrays.append((numpy.copy(array), cropval))  # Copy the array into a list
            del array
            band = ""
            img = ""
        i += 1

    count = 0
    finals = []

    for array, cropval in arrays:
        ltarrays = []
        nodataarrays = []
        for i in range(len(arrays)):
            if not i == count:
                # This bit finds the lowest value for each pixel, and something else I can't quite ascertain with the nodata
                lt = array.__lt__(arrays[i][0])
                ltarrays.append(numpy.copy(lt))
                ndarray = numpy.copy(array)
                ndarray[ndarray != nodata] = 0
                ndarray = numpy.rint(ndarray)
                ndarray = ndarray.astype(int)
                nodataarrays.append(ndarray)
        count += 1

        for i in range(len(ltarrays)):
            if i == 0:
                allpxbestfit = numpy.copy(ltarrays[i])
            else:
                allpxbestfit = allpxbestfit.__and__(ltarrays[i])
        finals.append(cropval * allpxbestfit)

    nodataarray = ""
    for ndarray in nodataarrays:
        if nodataarray == "":
            nodataarray = numpy.copy(ndarray)
        else:
            nodataarray = nodataarray.__and__(ndarray)  # Merge all nodata pixels from each fit image into one array

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
        searchdict["other"] = len(
            [x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
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

    searchdict["other"] = len(
        [x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
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

    return accuracy, classification, outstring


def classify_and_assess_accuracy(searchdir, cropimgpath, searchstringsvals, nodata,
                                 threshstart=500, threshstep=100, threshstepcount=10, outputdir=None,
                                 classifiedimagename=None, singlethresh=None,
                                 plotcorrectpx=False, plotincorrectpx=False):
    """
    """
    #TODO Docstring

    if outputdir is None:
        outputdir = searchdir
    else:
        pass

    if classifiedimagename is None:
        today = dt.now()
        classifiedimagename = today.strftime("%Y-%m-%d_%H%M_") + os.path.splitext(os.path.basename(cropimgpath))[0]

    classificationimage = os.path.join(outputdir, classifiedimagename + ".tif")
    accuracyreport = os.path.join(outputdir, classifiedimagename + ".txt")

    if plotcorrectpx:
        correctpxplot = Plot(outputdir, "correctpxplot")

    if plotincorrectpx:
        incorrectpxplot = Plot(outputdir, "incorrectpxplot")

    try:
        #np.set_printoptions(threshold=np.nan)  # For debug: Makes numpy print whole contents of an array.
        #Crop image is constant for all iterations
        cropimg = gdalObject()
        cropimg.open(cropimgpath)
        band = cropimg.gdal.GetRasterBand(1)
        croparray = band.ReadAsArray(0, 0, cropimg.cols, cropimg.rows)
        band = None
        cropimg.close()

        filelist = []
        files = os.listdir(searchdir)
        for f in files:
            for string, val in searchstringsvals:
                if f.endswith(".tif"):
                    if string in f:
                        filelist.append((os.path.join(searchdir, f), val))

        if singlethresh:
            thresholds = []
            for val in range(threshstart, (threshstepcount * threshstep + threshstart), threshstep):
                threshtemp = []
                for i in range(len(filelist)):
                    threshtemp.append(val)
                thresholds.append(threshtemp)
        else:
            thresholds = generate_thresholds(threshstart, threshstep, threshstepcount, len(filelist))

        writestring = ""
        bestacc = 0
        bestthresh = ""
        for thresh in thresholds:
            start = dt.now()
            accuracy, classification, outstring = classify_with_threshold(croparray,
                                                                                      filelist,
                                                                                      searchdir, searchstringsvals,
                                                                                      thresh, nodata)
            writestring = writestring + outstring

            if accuracy > bestacc:
                bestacc = accuracy
                bestthresh = thresh

            elapsed = dt.now() - start
            toprint = [thresh, "{}:{}".format(elapsed.seconds, elapsed.microseconds), accuracy, bestacc, bestthresh]
            width = (6 * len(filelist))
            print "Thresh: {: <{width}}   Time: {}   Acc: {: <14}   Best: {: <14} at {}\r".format(*toprint, width=width),

    except Exception as e:
        print e

    finally:

        with open(accuracyreport, 'w') as text:
            text.write("Classification using curves from {0}".format(searchdir))
            text.write("{0}\nBest:\n{1} {2}".format(writestring, bestthresh, bestacc))

        print "\n", bestthresh, bestacc

        accuracy, classification, outstring = classify_with_threshold(croparray, filelist, searchdir,
                                                                                  searchstringsvals, bestthresh,
                                                                                  nodata)

        driver = gdal.GetDriverByName("ENVI")
        driver.Register()

        outds = driver.Create(classificationimage, cropimg.cols, cropimg.rows, 1, GDT_Int16)
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

    #TODO docstring

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
        for i in range(bands):
            band = img.GetRasterBand(i + 1)
            v = int(band.ReadAsArray(location[1], location[0], 1, 1))
            print "\tBand {0}: {1}".format(i + 1, v)
            values.append(v)
            band = None
        refs.append(values)
    img = None

    return refs


def write_refs_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):
    print "Writing pixel curves to output file:"
    print referencevalues

    #TODO: Rewrite this and next function with append on open, and integrate together to reduce redundant code

    #output individual pixel curves to file
    with open(os.path.join(outdir, cropname + postfix + "_points.ref"), "w") as f:
        if comment:
            f.write("//" + comment + "\n\n")
        point = 1
        for points in referencevalues:
            f.write("\nPoint {0}:\n".format(point))
            point += 1
            imgnumber = 1
            for val in points:
                doy = band_number_to_doy(imgnumber, startdoy, doyinterval)
                f.write("{0} {1}\n".format(doy, val))
                imgnumber += 1


def write_mean_ref_to_txt(cropname, referencevalues, startdoy, doyinterval, outdir, comment="", postfix=""):

    #TODO docstring

    print "Writing mean reference curve to output file:"
    print referencevalues

    #output mean pixel values to file
    with open(os.path.join(outdir, cropname + postfix + "_mean.ref"), "w") as f:
        if comment:
            f.write("//" + comment + "\n\n")
        meanvals = get_mean_values(referencevalues)
        imgnumber = 1
        for val in meanvals:
            doy = band_number_to_doy(imgnumber, startdoy, doyinterval)
            f.write("{0} {1}\n".format(doy, val))
            imgnumber += 1


def get_mean_values(referencevalues):

    #TODO docstring

    array = numpy.array(referencevalues)
    mean = numpy.mean(array, axis=0)
    return list(mean)


def get_reference_curves(image, refstoget, startdoy, imageinterval, outdir="", filepostfix=""):
    """

    """

    #TODO docstring

    if not outdir:
        outdir = os.path.dirname(image)

    for shapefile in refstoget:
        cropname = os.path.splitext(os.path.basename(shapefile))[0]
        locs = get_px_coords_from_points(image, shapefile)
        refvals = get_crop_pixel_values(image, locs)
        comment = "Generated from {0} by get_crop_pixel_values version 0.1.".format(image)
        write_refs_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment, postfix=filepostfix)
        write_mean_ref_to_txt(cropname, refvals, startdoy, imageinterval, outdir, comment=comment,
                              postfix=filepostfix)


########## PROCEDURE ##########

if __name__ == '__main__':
    sys.exit()

