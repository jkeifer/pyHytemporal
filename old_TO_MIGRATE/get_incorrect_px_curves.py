import numpy as np
import os
from osgeo import gdal
from osgeo.gdalconst import *
import sys
from create_rule_image_multiprocessed_bypx import get_sort_dates_values
import build_multiband_image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import floor
from pyhytemporal.utils import create_output_dir


def color_picker(i):
    colorcodes = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = i % len(colorcodes)

    return colorcodes[i]

def plot_pixel(pixelproperties, dims=(11, 8), color=None):
    x, y = get_sort_dates_values(pixelproperties[4])
    figure = plt.figure()
    axes = figure.add_subplot(1,1,1)
    axes.plot(x, y, '{0}-'.format(color))
    axes.set_xlabel("Pixel col {0} row {1}: {2} as {3} in classified".format(pixelproperties[0], pixelproperties[1], pixelproperties[3], pixelproperties[2]))

    return figure

def get_incorrect_px_curves(inputpath, referencepath, vitimeseriespath, startdoy, interval, outputpath=None, nodata=None):

    if not outputpath:
        root = os.path.dirname(inputpath)
        outputpath = create_output_dir(root, "plots")

    #Open images
    inputimage = build_multiband_image.open_image(inputpath, returnimage=True)
    referenceimage = build_multiband_image.open_image(referencepath, returnimage=True)
    timeseries = build_multiband_image.open_image(vitimeseriespath, returnimage=True)

    #Get band 1 of input and reference images and read as arrays
    inband = inputimage[5].GetRasterBand(1)
    inarray = inband.ReadAsArray(0, 0, inputimage[0], inputimage[1])
    inband = ""

    refband = referenceimage[5].GetRasterBand(1)
    refarray = refband.ReadAsArray(0, 0, referenceimage[0], referenceimage[1])
    refband = ""

    #Close in and ref images

    #Find unique values in the input image
    uniquevals = set(np.unique(inarray).tolist())
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

    plotvalues = []
    i = 0
    for group, pixels in coordinatelist.items():
        for pixel in pixels:
            pxvalues = {}
            for key, val in arrays.items():
                doy = interval * (key - 1) + startdoy

                if doy >= 366:
                    t = floor((doy - 366)/interval)
                    doy = interval * t + 366

                pxvalues[doy] = val[pixel[1], pixel[0]]
            color = color_picker(i)
            pixel.append(pxvalues)
            pixel.append(color)
            plotvalues.append(pixel)
        i += 1

    pdf = PdfPages(os.path.join(outputpath, "plots.pdf"))
    for pixel in plotvalues:
        figure = plot_pixel(pixel, outputpath, color=pixel[5])
        pdf.savefig()
        plt.close(figure)
    pdf.close()



if __name__ == '__main__':

    inimage = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/test1_envicurves/fullpxonly/clip1refs/KansasNDVI_2012_clip1_SLSQP/classified_3crop.tif"
    refimage = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif"
    timeseries = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip1.tif"
    startdoy = 17
    interval = 16
    nodata = -3000

    sys.exit(get_incorrect_px_curves(inimage, refimage, timeseries, startdoy, interval, nodata=nodata))