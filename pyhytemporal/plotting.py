import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_color_picker(i):
    colorcodes = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = i % len(colorcodes)

    return colorcodes[i]

def plot_pixel(pixelObject, color=None):
    pdf = PdfPages(os.path.join(outputpath, "plots.pdf"))

    figure = plt.figure()
    axes = figure.add_subplot(1,1,1)
    axes.plot(pixelObject.values, bandDOYs, '{0}-'.format(color))
    axes.set_xlabel("Pixel row {0} col {1}: {2} as {3} in classified".format(pixelObject.row,
                                                                             pixelObject.col,
                                                                             pixelObject.actualclass,
                                                                             pixelObject.classifiedclass))

    pdf.savefig()
    plt.close(figure)
    pdf.close()


class pixel_plot(object):
    def __init__(self):
        self.plot = None
    def create_plot(self):

    def add_pixel_plot(pixelObject, color=None):

        pdf = PdfPages(os.path.join(outputpath, "plots.pdf"))

        figure = plt.figure()
        axes = figure.add_subplot(1,1,1)
        axes.plot(pixelObject.values, bandDOYs, '{0}-'.format(color))
        axes.set_xlabel("Pixel row {0} col {1}: {2} as {3} in classified".format(pixelObject.row,
                                                                                 pixelObject.col,
                                                                                 pixelObject.actualclass,
                                                                                 pixelObject.classifiedclass))

        pdf.savefig()
        plt.close(figure)
        pdf.close()

    def close_plot
