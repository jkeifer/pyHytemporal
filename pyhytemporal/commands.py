import os
import sys
import click
from utils import *
from classification import *


def validate_value(ctx, param, value):
    """
    Check to make sure the arg is fomatted correctly...
    """

    #TODO: Write this function

    return True

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
#TODO Add an option to use geographic or pixel extent (done) to clip raster in addition to clip to shape option
@click.option('-P', '--prompt-mode', is_flag=True, is_eager=True, expose_value=False, callback=find_fit_prompt,
              help="Enable prompt mode. This will prompt you for each of the arguments to the function. Use if you aren't good at the command line.")
def find_fit(vi, signaturedirectory, image, outputdir, outputfoldername, startdoy, doyinterval, temporalshift, threshold,
             ndvalue, subset, meantype, numberofprocesses, cliptopixelextent, cliptoshapeextent):
    """
    Fit the fit of reference temporal signatures to pixels in a multidate image.
    """
    #TODO Docstring
    #TODO Add Parameter Validation Callbacks as necessary

    signatures = signatureCollection(vi)
    sigFiles = find_files(signaturedirectory, "mean.ref", recursive=False)

    for f in sigFiles:
        signatures.add(f)

    if outputdir is None:
        outputdir = os.path.dirname(image)

    outdir = create_output_dir(outputdir, outputfoldername)

    if cliptoshapeextent:
        #TODO: Need to write a method to clip to the extent of a shapefile
        imagetoprocess = image
    elif cliptopixelextent:
        imagename, ext = os.path.splitext(os.path.basename(image))
        outimage = os.path.join(outdir, imagename + "_clip" + ext)
        imagetoprocess = clip_raster_to_extent(image, outimage, cliptopixelextent[0], cliptopixelextent[1],
                                               cliptopixelextent[2], cliptopixelextent[3])
    else:
        imagetoprocess = image

    phenological_classificaion(imagetoprocess, outdir, signatures, startdoy, doyinterval,
                               temporalshift, threshold=threshold, ndvalue=ndvalue, subset=subset, meantype=meantype,
                               workers=numberofprocesses)


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
              default=None, help="Path to the output directory. Default is to use the directory containing the image.")
@click.option('-s', '--startDOY', type=click.INT, help="The start DOY for the multidate image.", required=True)
@click.option('-d', '-DOYinterval', type=click.INT, help="The interval of the imagery in the multidate image.",
              required=True)
@click.option('-l', '--filelabel', type=click.STRING, default="",
              help="A label to postfix on each of the .ref file names")
def extract_signatures(image, shapefiledirectory, startdoy, doyinterval, outputdir, filelabel):
    """
    Extracts temporal signatures for a set of point geometry shapefiles in a specified directory and outputs them to a
    set of .ref files in an output directory.
    """

    if outputdir is None:
        outputdir = os.path.dirname(image)

    shapefiles = find_files(shapefiledirectory, ".shp", recursive=False)

    #TODO: Need a method to find only valid shapefiles in the directory

    get_reference_curves(image, shapefiles, startdoy, doyinterval, outdir=outputdir, filepostfix=filelabel)

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
@click.option('-n', '--ndvalue', type=click.INT, default=-3000,
              help="The value for NODATA in the multidate image and output fit images. Default is -3000.")
@click.option('-O', '--outputimagename', type=click.STRING, default=None,
              help="Name of the image to be created with the file extension. Default is the date and crop image name.")
def classify(fitimagedirectory, cropimage, outputdirectory, ndvalue, outputimagename, valueofcropinimage):
    """
    Classify a multidate image and assess the accuracy of said classification.
    """

    if outputdirectory is None:
        outputdirectory = os.path.dirname(fitimagedirectory)

    classify_and_assess_accuracy(fitimagedirectory, cropimage, valueofcropinimage, ndvalue,
                                 outdir=outputdirectory, outfilename=outputimagename)


cli.add_command(find_fit)
cli.add_command(build_multidate_image)
cli.add_command(extract_signatures)

if __name__ == '__main__':
    cli()