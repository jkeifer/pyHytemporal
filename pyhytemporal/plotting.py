import os
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 4}

rc('font', **font)


class FigureError(Exception):
    """
    Error for figure opening/closing in plotting lib
    """
    pass


class LegendError(Exception):
    """
    Error for figure opening/closing in plotting lib
    """
    pass


def plot_color_picker(i):
    from matplotlib.colors import cnames
    from random import shuffle
    colorcodes = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colorcodes.extend(cnames.keys())
    i = i % len(colorcodes)

    return colorcodes[i]


class Plot(object):
    def __init__(self, outputpath, namewithext):
        #TODO: Validate Directory
        self.pdf = PdfPages(os.path.join(outputpath, namewithext))
        self.figure = None
        self.lines = []

    def close_plot(self):
        try:
            self.close_figure()
        except FigureError:
            pass

        self.pdf.close()

    def create_figure(self):
        if not self.figure:
            self.figure = plt.figure()
        else:
            raise FigureError("Error: figure already exists. Close the figure to create a new figure.")

    def close_figure(self, save=True):
        if self.figure:
            if save:
                self.pdf.savefig()
            plt.close(self.figure)
            self.figure = None
        else:
            raise FigureError("Error: no figure is open. Cannot close nothing.")

    def add_legend(self):
        if self.lines:
            plt.legend()#labels=self.labels)
        else:
            raise LegendError("Error: No labels to add to the legend.")


class PixelPlot(Plot):
    def add_pixel(self, pixelObject, color=None, closefigure=False):

        if color:
            color = '{0}'.format(color)
        elif color is None and pixelObject.color:
            color = '{0}'.format(pixelObject.color)
        else:
            color = "black"

        try:
            self.create_figure()
        except FigureError:
            pass

        axes = self.figure.add_subplot(1, 1, 1)
        print("Plotting {0}, {1}".format(pixelObject.values, pixelObject.bandDOYs))
        axes.plot(pixelObject.bandDOYs, pixelObject.values, color=color, linestyle='-')
        axes.set_xlabel("Pixel row {0} col {1}: {2} as {3} in classified".format(pixelObject.row,
                                                                                 pixelObject.col,
                                                                                 pixelObject.actualclass,
                                                                                 pixelObject.classificationclass))
        if closefigure is True:
            self.close_figure()


class SignaturePlot(Plot):
    def plot_signature(self, temporalSignature, color=None, closefigure=False):

        if not color:
            color = "k"

        try:
            self.create_figure()
        except FigureError:
            pass

        axes = self.figure.add_subplot(1, 1, 1)
        self.lines.append(axes.plot(temporalSignature.daysofyear, temporalSignature.vivalues,
                                    color=color, linestyle='-', label=temporalSignature.name))
        #axes.set_xlabel("Temporal signature {0}".format(temporalSignature.name))

        if closefigure is True:
            self.close_figure()

    def plot_collection(self, signatureCollection, insamefigure=True, legend=True, closeplot=True):

        if insamefigure:
            closefig = False
        else:
            closefig = True

        for i, signature in enumerate(signatureCollection.signatures):
            color = plot_color_picker(i)
            self.plot_signature(signature, color=color, closefigure=closefig)

        if legend:
            self.add_legend()

        if closeplot:
            self.close_plot()


def plot_correct_incorrect_px(accuracyimage, multidateimageclipped, cropimage, classifiedimage, outputdir, plotcorrectpx=True, plotincorrectpx=True):
    # TODO: Finish below
    """Pesudocode:

    open accuracyimage as array
    open multidateimageclipped as array
    open cropimage as array
    open classifiedimage as array

    for row in accuracyimage
        for col in accuracyimage
            if accuracyimage[row, col] == 1 and plotcorrectpx:



    okay, need to create pixel objects for each non-nodata pixel in the accuracy image, get truth value and classified id,
    save into two lists, one of correct and one of incorrect
    plot each list"""

    if plotcorrectpx:
        correctpxplot = PixelPlot(outputdir, "correctpxplot")

    if plotincorrectpx:
        incorrectpxplot = PixelPlot(outputdir, "incorrectpxplot")