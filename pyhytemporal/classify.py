from datetime import datetime as dt
import os
import sys
import numpy
from osgeo import gdal
from osgeo.gdalconst import *
from core import gdalProperties
from imageFunctions import openImage, read_image_into_array, copySchemaToNewImage


OTHER_VALUE = 0

def classify_and_assess_accuracy(outputdir, cropimgpath, searchstringsvals, filevalist, nodata, thresholdlist,
                                 classifiedimagename=None, numberofprocesses=4):
    """
    """
    #TODO Docstring

    #from multiprocessing import Pool
    #pool = Pool(numberofprocesses)

    if classifiedimagename is None:
        today = dt.now()
        classifiedimagename = today.strftime("%Y-%m-%d_%H%M_") + os.path.splitext(os.path.basename(cropimgpath))[0]

    classificationimage = os.path.join(outputdir, classifiedimagename + ".tif")
    accuracyimage = os.path.join(outputdir, classifiedimagename + "_accuracy.tif")
    accuracyreport = os.path.join(outputdir, classifiedimagename + ".txt")

    #np.set_printoptions(threshold=np.nan)  # For debug: Makes numpy print whole contents of an array.
    #Crop image is constant for all iterations
    cropimg = openImage(cropimgpath)
    cropimgproperties = gdalProperties(cropimg)
    croparray = read_image_into_array(cropimg)
    cropimg = None

    arraylist = [(read_image_into_array(openImage(f[0])), f[1]) for f in filevalist]

    writestring = ""
    bestacc = 0
    bestthresh = None

    try:
        if len(thresholdlist) == 1:
            writestring = "\n\n**Only using a single threshold value--not iterating.**\n\n"
            bestthresh = thresholdlist[0]
        else:

            #TODO: Refactor to allow use of multiprocessing.Pool.map -- need to reason about the output/logging
            for thresh in thresholdlist:
                start = dt.now()
                accuracy, classification, outstring = classify_with_threshold(croparray, arraylist,
                                                                              searchstringsvals, thresh, nodata)
                writestring = writestring + outstring

                if accuracy > bestacc:
                    bestacc = accuracy
                    bestthresh = thresh

                elapsed = dt.now() - start
                toprint = [thresh, "{}:{}".format(elapsed.seconds, str(elapsed.microseconds).zfill(6)), accuracy, bestacc, bestthresh]
                width = (6 * len(arraylist))
                sys.stdout.write("Thresh: {: <{width}}   Time: {}   Acc: {: <14}   Best: {: <14} at {}\r".format(*toprint, width=width))
                sys.stdout.flush()

    except Exception as e:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print e
        traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stdout)

    finally:

        accuracy, classificationarray, outstring = classify_with_threshold(croparray, arraylist,
                                                                      searchstringsvals, bestthresh, nodata)

        writestring = writestring + outstring

        accuracyarray = find_correct_incorrect_array(croparray, classificationarray, ndvalue=nodata)

        # TODO: This repeat smells. Fix it.
        with open(accuracyreport, 'w') as text:
            text.write("Classification using fit images from {0}".format(os.path.dirname(filevalist[0][0])))
            text.write("{0}\nBest:\n{1} {2}".format(writestring, bestthresh, accuracy))

        print("\n{0}, {1}".format(bestthresh, accuracy))

        driver = gdal.GetDriverByName("ENVI")
        driver.Register()

        write_output_image(cropimgproperties, classificationimage, classificationarray, nodata)
        write_output_image(cropimgproperties, accuracyimage, accuracyarray, nodata)

        print("outputted")

        return 0


def classify_with_threshold(croparray, arraylist, searchstringsvals, thresh, nodata):
    """
    """
    #TODO DOCSTRING
    #TODO test to ensure thresh length is equal to the number of image files?
    #TODO REFACTOR INTO SMALLER, MORE-TESTABLE FUNCTIONS!!!!

    arrays = []

    for i, arrayval in enumerate(arraylist):
        array = numpy.copy(arrayval[0])  # Not sure if I need to cpy, but believe the array is passed to the func by ref
        #TODO: Test to see if arrays are passed by reference
        array[array > thresh[i]] = 10000  # Make all values in the array greater than the thresh equal 10000
        arrays.append((array, arrayval[1]))  # Copy the array into a list

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


def find_correct_incorrect_array(trutharray, classificationarray, ndvalue=-3000):
    uniquevals = set(classificationarray)

    if OTHER_VALUE in uniquevals:
        uniquevals.remove(OTHER_VALUE)

    accuracyarray = numpy.array(classificationarray.__eq__(trutharray), dtype=numpy.int_)  # sets acc array to 1 where equal, 0 where not
    accuracyarray[classificationarray == ndvalue] = ndvalue  # where classification contains nodata, set acc to nodata
    accuracyarray[trutharray == ndvalue] = ndvalue  # where trutharray contains nodata, set acc to nodata
    accuracyarray[classificationarray == OTHER_VALUE and not trutharray in uniquevals] = 1  # where classified is other and truth is other, set to correct (1)
    return accuracyarray


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


def get_fit_rasters(searchdir, searchstringsvals):
    #Find fit images and open as arrays, building a list of tuples of the array and crop value
    files = os.listdir(searchdir)
    filevallist = [(os.path.join(searchdir, f), val) for string, val in searchstringsvals
                for f in files if f.endswith(".tif") and string in f]
    return filevallist


def write_output_image(propertiestocopy, outname, data, nodata):
    raster = copySchemaToNewImage(propertiestocopy, outname, datatype=GDT_Int16)
    raster_band = raster.GetRasterBand(1)
    raster_band.WriteArray(data, 0, 0)
    raster_band.SetNoDataValue(nodata)
    raster_band.FlushCache()
    raster_band = None
    raster = None