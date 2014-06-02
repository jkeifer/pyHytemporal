import os
import sys
import click
from utils import *
from classification import *

def find_fit_prompt(ctx, value):
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
@click.option('-n', '--numberofprocesses', type=click.INT, default=4,
              help="The max number of processes to spawn at any time. Default is 4. Set lower or higher depending on number of processors/cores in your machine.")
@click.option('-c', '--cliptoshapeextent', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,
                                                resolve_path=True), default=None,
              help="Path to a shapefile to clip the raster. If omitted, entire raster extent will be processed.")
@click.option('-C', '--cliptopixelextent', nargs=4, type=int, default=None,
              help="Pixel coordinates and number of pixels to clip raster. For example, entering \"2482, 1089, 100, 100\"\
                    will create a 100px square image going right and down from pixel 2482, 1089 in the original image.")
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

if __name__ == '__main__':
    find_fit()