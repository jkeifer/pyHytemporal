__author__ = 'phoetrymaster'

import create_rule_image_multiprocessed_bypx
import matplotlib.pyplot as plt
import matplotlib.legend as legend
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

def main():
    outpath = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/GIS In Action/Images/example_fit_process_1.pdf"
    wheat = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012clip1test2/clip1/wwheat_NDVI_mean.ref"
    corn = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012clip1test2/clip1/corn_NDVI_mean.ref"
    soy = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012clip1test2/clip1/soy_NDVI_mean.ref"

    pixel = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Scripting/Testing/corn_example_point.ref"

    wheatmean = create_rule_image_multiprocessed_bypx.read_reference_file(wheat)[1]
    cornmean = create_rule_image_multiprocessed_bypx.read_reference_file(corn)[1]
    soymean = create_rule_image_multiprocessed_bypx.read_reference_file(soy)[1]

    refs = [('wheat', wheatmean), ('corn', cornmean), ('soy', soymean)]
    print refs

    pixelvalues = create_rule_image_multiprocessed_bypx.read_reference_file(pixel)[1]

    toplot = [(wheatmean, 'b', 'Wheat Signature'), (cornmean, 'y', "Corn Signature"), (soymean, 'g', 'Soy Signature'), (pixelvalues, 'r', 'Unknown Pixel')]

    figure = plt.figure()
    axes = figure.add_subplot(1, 1, 1)

    pdf = PdfPages(outpath)

    for item, color, label in toplot:
        print item
        x, y = create_rule_image_multiprocessed_bypx.get_sort_dates_values(item)

        newy = []
        for value in y:
            newy.append(value / 10000.0)

        axes.plot(x, newy, '{0}-'.format(color), label=label)

    handles, labels = axes.get_legend_handles_labels()

    axes.set_xlabel("Day of Year")
    axes.set_ylabel("NDVI")
    axes.legend(handles, labels, fontsize='small')
    pdf.savefig()
    plt.close(figure)
    pdf.close()

    for crop, ref in refs:
        interp = create_rule_image_multiprocessed_bypx.interpolate_ref(ref)
        result = create_rule_image_multiprocessed_bypx.find_fit(pixelvalues, interp, 0, 0, 'SLSQP', bnds=((0.6, 1.4), (0.6, 1.4), (-10, 10)))
        print(crop + ":")
        print(result)


if __name__ == '__main__':
    sys.exit(main())