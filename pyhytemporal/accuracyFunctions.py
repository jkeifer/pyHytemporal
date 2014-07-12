from imageFunctions import openImage
import numpy

def get_incorrect_pixels(inputpath, referencepath, vitimeseriespath, startdoy, interval, nodata=None):

    #Open images
    inputimage = openImage(inputpath)
    referenceimage = openImage(referencepath)
    timeseries = openImage(vitimeseriespath)

    #Get band 1 of input and reference images and read as arrays
    inband = inputimage[5].GetRasterBand(1)
    inarray = inband.ReadAsArray(0, 0, inputimage[0], inputimage[1])
    inband = ""

    refband = referenceimage[5].GetRasterBand(1)
    refarray = refband.ReadAsArray(0, 0, referenceimage[0], referenceimage[1])
    refband = ""

    #Close in and ref images

    #Find unique values in the input image
    uniquevals = set(numpy.unique(inarray).tolist())
    uniquevals.remove(nodata)
    print uniquevals

    #Iterate through pixels in images, comparing arrays for each value in input image, making note of tested values and excluding them from comparison with 0 value in input image
    #Return pixel coordinates and the identified and true identities of the pixels from the input and reference images
    coordinatelist = {}
    for value in uniquevals:
        coordinatelist[value] = []

    numincorrect = 0

    for col in range(0, inputimage[0]):
        for row in range(0, inputimage[1]):

            if (inarray[row, col] != 0) and (inarray[row, col] != refarray[row, col]) and (inarray[row, col] != nodata):

                if refarray[row, col] in uniquevals:
                    coordinatelist[refarray[row, col]].append([col, row, inarray[row, col], refarray[row, col]])
                    numincorrect += 1
                else:
                    coordinatelist[0].append([col, row, inarray[row, col], refarray[row, col]])
                    numincorrect += 1

            elif (inarray[row,col] == 0) and (refarray[row, col] in uniquevals):
                coordinatelist[refarray[row, col]].append([col, row, inarray[row, col], refarray[row, col]])
                numincorrect += 1

    #For each returned pixel, extract and plot the curve from vitimeseries and export plot as image; be sure to label plot correctly
    bandcnt = timeseries[5].RasterCount
    arrays = {}
    for i in range(1, bandcnt + 1):
        band = timeseries[5].GetRasterBand(i)
        data = band.ReadAsArray(0, 0, timeseries[0], timeseries[1])
        arrays[i] = data

    print("Found {0} incorrect pixels to process...".format(numincorrect))
