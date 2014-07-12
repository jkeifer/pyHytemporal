import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#from core import signatureCollection, temporalSignature


def plot_color_picker(i):
    colorcodes = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = i % len(colorcodes)

    return colorcodes[i]


class Plot(object):
    def __init__(self, outputpath, name):
        self.pdf = PdfPages(os.path.join(outputpath, name + ".pdf"))
        self.figure = None

    def create_figure(self):
        if not self.figure:
            self.figure = plt.figure()
        else:
            raise Exception("Error: figure already exists. Close the figure to create a new figure.")

    def close_figure(self, save=True):
        if self.figure:
            if save:
                self.pdf.savefig()
            plt.close(self.figure)
        else:
            raise Exception("Error: no figure is open. Cannot close nothing.")

    def close_plot(self):
        self.pdf.close()


class PixelPlot(Plot):
    def add_pixel_plot(self, pixelObject, color=None, closefigure=False):

        if color:
            color = '{0}-'.format(color)
        elif color is None and pixelObject.color:
            color = '{0}-'.format(pixelObject.color)
        else:
            color = "black-"

        self.create_figure()
        axes = self.figure.add_subplot(1, 1, 1)
        axes.plot(pixelObject.values, pixelObject.bandDOYs, color)
        axes.set_xlabel("Pixel row {0} col {1}: {2} as {3} in classified".format(pixelObject.row,
                                                                                 pixelObject.col,
                                                                                 pixelObject.actualclass,
                                                                                 pixelObject.classifiedclass))
        if closefigure is True:
            self.close_figure()


class SignaturePlot(Plot):
    def plot_signature(self, temporalSignature, color=None, closefigure=False):

        if color:
            color = '{0}-'.format(color)
        else:
            color = "black-"

        self.create_figure()
        axes = self.figure.add_subplot(1,1,1)
        axes.plot(temporalSignature.vivalues, temporalSignature.daysofyear, color)
        axes.set_xlabel("Temporal signature {1}".format(temporalSignature.name))

        if closefigure is True:
            self.close_figure()

    def plot_collection(self, signatureCollection, insamefigure=True):

        if insamefigure:
            closefig = False

        for i, signature in enumerate(signatureCollection.signatures):
            color = plot_color_picker(i)
            self.plot_signature(signature, color=color, closefigure=closefig)

