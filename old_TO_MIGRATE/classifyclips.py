from create_rule_image_multiprocessed_bypx import phenological_classificaion, read_reference_file
import os
from pyhytemporal.utils import find_files

clip1refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip1", "mean.ref")
clip2refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip2", "mean.ref")
clip3refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip3", "mean.ref")
clip4refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip4", "mean.ref")
clip5refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip5", "mean.ref")
clip6refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012/clip6", "mean.ref")
meanrefs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012", ".ref", recursive=False)

reffiles = [clip1refs, clip2refs, clip3refs, clip4refs, clip5refs, clip6refs, meanrefs]
rootout = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified"
outfolders = [os.path.join(rootout, "clip1refs"), os.path.join(rootout, "clip2refs"), os.path.join(rootout, "clip3refs"), os.path.join(rootout, "clip4refs"), os.path.join(rootout, "clip5refs"), os.path.join(rootout, "clip6refs"), os.path.join(rootout, "meanrefs")]

clip1imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip1.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"]
clip2imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip2.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip2.tif"]
clip3imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip3.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip3.tif"]
clip4imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip4.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip4.tif"]
clip5imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip5.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip5.tif"]
clip6imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip6.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip6.tif"]

imagelist = [clip1imgs, clip2imgs, clip3imgs, clip4imgs, clip5imgs, clip6imgs]

searchstrings = ["soy", "corn", "wwheat", "sorghum", "wwheatsoydbl"]
fitmethods = ["SLSQP"]#, "TNC"]


for method in fitmethods:
    for images in imagelist:
        for img in images:
            name = os.path.splitext(os.path.basename(img))[0]
            if "NDVI" in name:
                type = "NDVI"
            else:
                type = "EVI"

            for reffile, outfolder in zip(reffiles, outfolders):

                refs = {}
                for f in reffile:
                    if type in f:
                        for string in searchstrings:
                            if string + "_" + type in f:
                                error, refs[string] = read_reference_file(f)
                if not error:
                    print refs
                    phenological_classificaion(img, outfolder, name + "_" + method, refs, "ENVI", 17, 16, -500, 0, method, toprint=False)
                else:
                    print "ERROR +++++++++++++++++++++++++++++++++++++++++++++++++++ ERROR"