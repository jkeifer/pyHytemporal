import numpy
import os
from osgeo import gdal
from osgeo.gdalconst import *
import sys
from datetime import datetime as dt

from itertools import product

def generate_thresholds(start, step, numberofsteps, lengthofelement):
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

    u = 0
    for f, cropval in filelist:
        img = gdal.Open(f, GA_ReadOnly)
        if img is None:
            raise Exception("Could not open: {0}".format(os.path.join(searchdir, f)))
        else:
            rows = img.RasterYSize
            cols = img.RasterXSize
            band = img.GetRasterBand(1)
            array = band.ReadAsArray(0, 0, cols, rows)
            array[array > thresh[u]] = 10000
            #print thresh[u]
            arrays.append((numpy.copy(array), cropval))
            band = ""
            img = ""
        u += 1

    count = 0
    finals = []

    #with open(os.path.join(searchdir, "arrays.txt"), 'w') as text:
    #    for i, j in arrays:
    #        text.write(str(i)+"\n\n")

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

        #with open(os.path.join(searchdir, "ltarrays.txt"), 'w') as text:
        #    for i in ltarrays:
        #        text.write(str(i)+"\n\n")

        for i in range(0, len(ltarrays)):
            if not i:
                allpxbestfit = numpy.copy(ltarrays[i])
            else:
                allpxbestfit = allpxbestfit.__and__(ltarrays[i])
        finals.append(cropval * allpxbestfit)

    #with open(os.path.join(searchdir, "finals.txt"), 'w') as text:
    #    for i in finals:
    #        text.write(str(i))

    nodataarray = ""
    for ndarray in nodataarrays:
        if nodataarray == "":
            nodataarray = numpy.copy(ndarray)
        else:
            nodataarray = nodataarray.__and__(ndarray)

    #print nodataarray

    classification = ""
    for final in finals:
        if classification == "":
            classification = numpy.copy(final)
        else:
            classification = classification.__or__(final)
    classification = classification.__or__(nodataarray)
    #print classification

    #with open(os.path.join(searchdir, "classified.txt"), 'w') as text:
    #    text.write(str(classification))

    #Accuracy Assessment
    results = {}
    for string, val in searchstringsvals:
        dict = {}
        temparray = numpy.copy(classification)
        temparray[temparray != val] = 0
        correct = temparray.__eq__(croparray)
        incorrect = temparray.__ne__(croparray)
        temparray[temparray == val] = 1
        incorrectvals = incorrect.__mul__(croparray).__mul__(temparray)
        for string2, val2 in searchstringsvals:
            temparray2 = numpy.copy(incorrectvals)
            temparray2[temparray2 != val2] = 0
            dict[string2] = temparray2.sum() / val2

        dict[string] = correct.sum()
        dict["other"] = len([x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
        results[string] = dict.copy()

    dict = {}
    temparray = numpy.copy(classification)
    temparray[classification == 0] = 1
    temparray[classification != 0] = 0
    incorrectvals = temparray.__mul__(croparray)

    for string2, val2 in searchstringsvals:
        temparray2 = numpy.copy(incorrectvals)
        temparray2[temparray2 != val2] = 0
        dict[string2] = temparray2.sum() / val2

    dict["other"] = len([x for y in incorrectvals for x in y if not x in zip(*searchstringsvals)[1] and not x == 0])
    results["other"] = dict.copy()

    numpx = 0
    for key, val in results.items():
        for k, v in val.items():
            numpx += v

    #print "\n"
    #print results
    #print numpx


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
        h = h + total
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
        cropimg = gdal.Open(cropimgpath, GA_ReadOnly)
        if cropimg is None:
            raise Exception("Could not open: {0}".format(cropimgpath))
        else:
            rows = cropimg.RasterYSize
            cols = cropimg.RasterXSize
            projection = cropimg.GetProjection()
            transformation = cropimg.GetGeoTransform()
            band = cropimg.GetRasterBand(1)
            datatype = band.DataType
            croparray = band.ReadAsArray(0, 0, cols, rows)
            band =""
            cropimg = ""
            print "Opened crop img"

        #
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

        accuracy, classification, cols, rows, outstring = classify_with_threshold(croparray,
                                                                                           filelist,
                                                                                           searchdir, searchstringsvals,
                                                                                           bestthresh, nodata)
        driver = gdal.GetDriverByName("ENVI")
        driver.Register()

        outds = driver.Create(outFile, cols, rows, 1, GDT_Int16)
        outds.SetGeoTransform(transformation)
        outds.SetProjection(projection)
        outband = outds.GetRasterBand(1)
        outband.WriteArray(classification, 0, 0)
        outband.SetNoDataValue(-3000)
        outband.FlushCache()

        outband = ""
        outds = ""

        print "outputted"

        return 0

if __name__ == '__main__':

    searchdir = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/test1_envicurves/fullpxonly/clip1refs/KansasNDVI_2012_clip1_SLSQP/"
    cropimgpath = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
    searchstringsvals = [("soy.", 5), ("wwheat.", 24), ("corn", 1)]#, ("sorghum", 4), ("wwheatsoydbl", 26)]
    nodata = -3000

    sys.exit(main(searchdir, cropimgpath, searchstringsvals, nodata))
