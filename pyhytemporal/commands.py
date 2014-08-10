import click


def validate_value(ctx, param, value):
    """
    Check to make sure the arg is fomatted correctly...
    """

    #TODO: Write this function
    toreturn = []
    for v in value:
        toreturn.append((str(v[0]), int(v[1])))

    return toreturn

@click.group()
def cli():
    pass


def find_fit_prompt(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    ctx.abort()


@click.command()
@click.option('-i', '--image', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                               resolve_path=True),
              required=True, help="Path to the multidate image file.")
@click.option('-j', '--signaturedirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True,
                                                            resolve_path=True),
              required=True, help="Path to the directory containing the temporal signature files to be used.")
@click.option('--vi', type=click.STRING, default="", help="Name of the VI being used. Default is a blank string.")
@click.option('-o', '--outputdir', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              default=None, help="Path to the output directory. Default is to use the directory containing the image.",)
@click.option('-f', '--outputfoldername', type=click.STRING, default='fit_images',
              help="Name of the folder to be created with the output files. Default is 'fit_images'.")
@click.option('-s', '--startDOY', type=click.INT, help="The start DOY for the multidate image.", required=True)
@click.option('-d', '-DOYinterval', type=click.INT, help="The interval of the imagery in the multidate image.",
              required=True)
@click.option('-t', '--temporalshift', type=click.INT, default=0,
              help="A guess to the temporal shift of the multidate image from the temporal signatures. Default is 0 days.")
@click.option('-T', '--threshold', type=click.INT, default=None,
              help="A value beneath which all VI values will be ignored. Default is none.")
@click.option('-n', '--ndvalue', type=click.INT, default=-3000,
              help="The value for NODATA in the multidate image and output fit images. Default is -3000.")
@click.option('-S', '--subset', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                                resolve_path=True),
              default=None, help="A shapefile of points for each pixel to be fit. Used to eliminate mixels. Default is none.")
@click.option('-m', '--meantype', type=click.Choice(['arithmetic', 'geometric']), default='arithmetic',
              help="The type of mean (arithmetic or geometric) used in the RMSE fitting of the signatures to the pixels. Default is arithmetic.")
@click.option('-p', '--numberofprocesses', type=click.INT, default=4,
              help="The max number of processes to spawn at any time. Default is 4. Set lower or higher depending on number of processors/cores in your machine.")
@click.option('-c', '--cliptoshapeextent', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                                resolve_path=True), default=None,
              help="Path to a shapefile to clip the raster. If omitted, entire raster extent will be processed.")
@click.option('-C', '--cliptopixelextent', nargs=4, type=int, default=None,
              help="Pixel coordinates and number of pixels to clip raster. For example, entering \"2482, 1089, 100, 100\" will create a 100px square image going right and down from pixel 2482, 1089 in the original image.")
@click.option('--timebounds', nargs=2, type=int, default=None,
              help="Number of days to allow curve shifting before and after: -10, 10 is default and allows the curve to be shifted 10 days in either direction.")
@click.option('--xbounds', nargs=2, type=float, default=None,
              help="Bonds of allowable x-scaling: default is 0.6 and 1.4, allowing the curve to be stretched horizontally between 60% and 140% of initial width.")
@click.option('--ybounds', nargs=2, type=float, default=None,
              help="Bonds of allowable y-scaling: default is 0.6 and 1.4, allowing the curve to be stretched vertically between 60% and 140% of initial height.")
#TODO Add an option to use geographic or pixel extent (done) to clip raster in addition to clip to shape option
@click.option('-P', '--prompt-mode', is_flag=True, is_eager=True, expose_value=False, callback=find_fit_prompt,
              help="**CURRENTLY DISABLED** Enable prompt mode. This will prompt you for each of the arguments to the function. Use if you aren't good at the command line.")
def find_fit(vi, signaturedirectory, image, outputdir, outputfoldername, startdoy, doyinterval, temporalshift,
             threshold, ndvalue, subset, meantype, numberofprocesses, cliptopixelextent, cliptoshapeextent, timebounds,
             xbounds, ybounds):
    """
    Fit the fit of reference temporal signatures to pixels in a multidate image.
    """
    #TODO Docstring
    #TODO Add Parameter Validation Callbacks as necessary

    # validate clip options
    if cliptopixelextent and cliptoshapeextent:
        click.BadParameter("Cannot clip the image to both a shapefile and pixel extent. Choose one or the other.")

    # import required modules
    import os
    from signatureFunctions import get_sigs_in_dir
    from utils import create_output_dir
    from imageFunctions import clip_raster_to_extent, clip_and_mask_raster_with_shapefile
    from fitting import fit_refs_to_image

    signatures = get_sigs_in_dir(signaturedirectory, viname=vi)

    if outputdir is None:
        outputdir = os.path.dirname(image)

    outdir = create_output_dir(outputdir, outputfoldername)

    if cliptoshapeextent:
        imagename, ext = os.path.splitext(os.path.basename(image))
        outimage = os.path.join(outdir, imagename + "_clip" + ext)
        imagetoprocess = clip_and_mask_raster_with_shapefile(image, cliptoshapeextent, outimage)
    elif cliptopixelextent:
        imagename, ext = os.path.splitext(os.path.basename(image))
        outimage = os.path.join(outdir, imagename + "_clip" + ext)
        imagetoprocess = clip_raster_to_extent(image, outimage, cliptopixelextent[0], cliptopixelextent[1],
                                               cliptopixelextent[2], cliptopixelextent[3])
    else:
        imagetoprocess = image

    fit_refs_to_image(imagetoprocess, outdir, signatures, startdoy, doyinterval,
                               temporalshift, threshold=threshold, ndvalue=ndvalue, subset=subset, meantype=meantype,
                               workers=numberofprocesses, timebounds=timebounds, xbounds=xbounds, ybounds=ybounds)


@click.command()
@click.option('-d', '--imagedirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True,
                                                            resolve_path=True),
              required=True, help="Path to the directory containing the .hdf image files to be used.")
@click.option('-n', '--outputimagename', type=click.STRING, default='multidate_image.tif',
              help="Name of the image to be created with the file extension. Default is 'multidate_image.tif'.")
@click.option('--vi', type=click.STRING, default="NDVI", help="Name of the VI to be used. Default is NDVI.")
@click.option('-o', '--outputdir', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              default=None, help="Path to the output directory. Default is to use the directory containing the image.",)
@click.option('-f', '--outputfoldername', type=click.STRING, default='multidate_image',
              help="Name of the folder to be created for the output file. Default is 'multidate_image'.")
@click.option('-N', '--ndvalue', type=click.INT, default=-3000,
              help="The value for NODATA in the multidate image and output fit images. Default is -3000.")
@click.option('-D', '--drivercode', type=click.STRING, default='GTiff',
              help="GDAL driver code for output image format. Default is GeoTIFF. Ensure output name extension is correct if using a different format.")
def build_multidate_image(imagedirectory, outputimagename, outputdir, outputfoldername, vi, drivercode, ndvalue):
    """
    Search directory for HDF MODIS files, get a VI from each HDF, and build single-date VI images in to a multi-date
    composite image.
    """
    from imageFunctions import build_multiband_image

    build_multiband_image(imagedirectory, outputimagename, outputfoldername, vi, str(drivercode), ndvalue,
                          outputdir=outputdir)


@click.command()
@click.option('-i', '--image', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                               resolve_path=True),
              required=True, help="Path to the multidate image file.")
@click.option('-v', '--shapefiledirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True,
                                                            resolve_path=True),
              required=True, help="Path to the directory containing point .shp files for each of the classes.")
@click.option('-o', '--outputdir', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              default=None, help="Path to the output directory. Default is to create a directory in the folder containing the image.")
@click.option('-s', '--startDOY', type=click.INT, help="The start DOY for the multidate image.", required=True)
@click.option('-d', '--DOYinterval', type=click.INT, help="The interval of the imagery in the multidate image.",
              required=True)
@click.option('-l', '--filelabel', type=click.STRING, default="",
              help="A label to postfix on each of the .ref file names")
@click.option('-l', '--filelabel', type=click.STRING, default="",
              help="A label to postfix on each of the .ref file names")
@click.option('-p', '--plotsigs', is_flag=True,
              help="Create a pdf plot of all the generated signatures.")
def extract_signatures(image, shapefiledirectory, startdoy, doyinterval, outputdir, filelabel, plotsigs):
    """
    Extracts temporal signatures for a set of point geometry shapefiles in a specified directory and outputs them to a
    set of .ref files in an output directory.
    """
    import os
    from plotting import SignaturePlot
    from utils import find_files, create_output_dir
    from signatureFunctions import get_sigs_in_dir, get_reference_curves

    if outputdir is None:
        outputdir = create_output_dir(os.path.dirname(image), "signatures", usetime=True)

    shapefiles = find_files(shapefiledirectory, ".shp", recursive=False)

    #TODO: Need a method to find only valid shapefiles in the directory

    get_reference_curves(image, shapefiles, startdoy, doyinterval, outdir=outputdir, filepostfix=filelabel)

    if plotsigs:
        sigs = get_sigs_in_dir(outputdir)
        plot = SignaturePlot(outputdir, "signaturePlot")
        plot.plot_collection(sigs)


@click.command()
@click.option('-i', '--fitimagedirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              required=True, help="Path to the directory containing the crop fit images.")
@click.option('-c', '--cropimage', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                               resolve_path=True),
              required=True, help="Path to the crop ground truth image file.")
@click.option('-o', '--outputdirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              default=None, help="Path to the output directory. Default is to use the directory containing the image.")
@click.option('-v', '--valueofcropinimage', multiple=True, nargs=2, callback=validate_value,
              help="The class name and its value in the crop image used for the accuracy assessment. E.g. \"Corn 1\"")
@click.option('-t', '--thresholds', default=[], type=click.STRING,
              help="A list of threshold values to use. Format each entry as a tuple in a python list with no spaces e.g. [(800,500,1200)]. Cannot be used with threshold stepping.")
@click.option('-n', '--ndvalue', type=click.INT, default=-3000,
              help="The value for NODATA in the multidate image and output fit images. Default is -3000.")
@click.option('-O', '--outputimagename', type=click.STRING, default=None,
              help="Name of the image to be created with the file extension. Default is the date and crop image name.")
@click.option('--tstart', type=click.INT,
              help="The threshold start value.")
@click.option('--tstep', type=click.INT,
              help="The threshold step value.")
@click.option('--tstepcount', type=click.INT,
              help="The number of threshold steps.")
@click.option('--nocombo', is_flag=True,
              help="Does not find combination of threshold steps, but steps through a single threshold value applied to all fit images.")
def classify(fitimagedirectory, cropimage, outputdirectory, ndvalue, outputimagename, valueofcropinimage, tstart, tstep,
             tstepcount, nocombo, thresholds):
    """
    Classify a multidate image and assess the accuracy of said classification.
    """

    # import required functions
    import os
    from utils import create_output_dir
    from classify import classify_and_assess_accuracy, generate_thresholds, get_fit_rasters

    # get the fit rasters to use
    filevallist = get_fit_rasters(fitimagedirectory, valueofcropinimage)

    # validate threshold parameters
    if (tstart or tstep or tstepcount or nocombo) and thresholds:
        raise click.BadParameter("Cannot use both a threshold list and stepping threshold options.")
    elif thresholds:
        thresholds = eval(thresholds)
        for thresh in thresholds:
            if len(thresh) != len(filevallist):
                raise click.BadParameter("Length of threshold in threshold value list is not the same as the number of fit rasters. Counts must be equal.")
            else:
                pass
        thresholds = (thresholds, len(thresholds))
    elif tstart and tstepcount and tstep:
        # create threshold generator
        if nocombo:
            thresholds = []
            for val in range(tstart, (tstepcount * tstep + tstart), tstep):
                threshtemp = [val for item in filevallist]
                thresholds.append(threshtemp)
            thresholds = (thresholds, len(thresholds))
        else:
            thresholds = (generate_thresholds(tstart, tstep, tstepcount, len(filevallist)), tstepcount**len(filevallist))
    else:
        raise click.BadParameter("Threshold options incomplete or otherwise incorrectly used.")

    if outputdirectory is None:
        outputdirectory = create_output_dir(os.path.dirname(fitimagedirectory), "classification", usetime=True)

    classify_and_assess_accuracy(outputdirectory, cropimage, valueofcropinimage, filevallist, ndvalue, thresholds,
                                 classifiedimagename=outputimagename)


@click.command()
@click.option('-i', '--multidateraster', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                               resolve_path=True),
              required=True, help="Path to the multidate raster file.")
@click.option('-p', '--pointfile', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                                resolve_path=True), default=None, required=True,
              help="Path to a point shapefile containing the points to be plotted")
@click.option('-s', '--startDOY', type=click.INT, help="The start DOY for the multidate image.", required=True)
@click.option('-d', '--DOYinterval', type=click.INT, help="The interval of the imagery in the multidate image.",
              required=True)
def plot_points(multidateraster, pointfile, startdoy, doyinterval):
    """

    """
    import os
    from utils import unique_name
    from plotting import PixelPlot
    from core import pixel as pixelObject
    from vectorFunctions import get_px_coords_from_shapefile
    from imageFunctions import openImage

    outpath = unique_name(os.path.dirname(multidateraster), "plots", ext=".pdf", usetime=True)

    coords = get_px_coords_from_shapefile(multidateraster, pointfile)

    plot = PixelPlot(os.path.dirname(outpath), os.path.basename(outpath))
    raster = openImage(multidateraster)

    for coord in coords:
        pixel = pixelObject(coord[0], coord[1])
        pixel.get_pixel_values(raster, startdoy, doyinterval)
        plot.add_pixel(pixel, closefigure=True)

    plot.close_plot()
    raster = None


@click.command()
@click.option('-j', '--signaturedirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True,
                                                            resolve_path=True),
              required=True, help="Path to the directory containing the temporal signature files to be used.")
@click.option('-o', '--outputdirectory', type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True,
                                                   readable=True, resolve_path=True),
              default=None, help="Path to the output directory. Default is to use the directory containing the signatures.")
@click.option('-n', '--name', type=click.STRING, default='signatures.pdf',
              help="Name of the plot pdf to be created with the file extension. Default is 'signatures.pdf'.")
def plot_sigs(signaturedirectory, outputdirectory, name):
    """

    """
    import os
    from utils import find_files, unique_name
    from core import signatureCollection
    from plotting import SignaturePlot

    if not outputdirectory:
        outputdirectory = signaturedirectory

    sigs = find_files(signaturedirectory, "mean.ref")

    if not sigs:
        click.BadParameter("Did not find any signature files in the specified directory.")

    signatures = signatureCollection()

    for sig in sigs:
        try:
            signatures.add(sig)
        except Exception as e:
            print e

        #TODO Fix core temporalSignature to use exceptions so they can be properly handled here

    name, ext = os.path.splitext(name)
    path = unique_name(outputdirectory, name, ext=ext)

    print("Outputting to {0}".format(path))

    plot = SignaturePlot(outputdirectory, os.path.basename(path))
    plot.plot_collection(signatures)


cli.add_command(find_fit)
cli.add_command(build_multidate_image)
cli.add_command(extract_signatures)
cli.add_command(classify)
cli.add_command(plot_points)
cli.add_command(plot_sigs)

if __name__ == '__main__':
    cli()