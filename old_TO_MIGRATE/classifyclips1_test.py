from create_rule_image_multiprocessed_bypx import phenological_classificaion, read_reference_file
import os
from get_px_coords_from_point import get_px_coords_from_points
from pyhytemporal.utils import find_files

reffiles = []
clip1refs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012clip1test2/clip1", "mean.ref")
reffiles.append(clip1refs)

#meanrefs = find_files("/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012", ".ref", recursive=False)
#reffiles.append(meanrefs)

rootout = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/test1_envicurves/fullpxonly/geometricmean"
#outfolders = [os.path.join(rootout, "clip1refs")]#, os.path.join(rootout, "clip2refs"), os.path.join(rootout, "clip3refs"), os.path.join(rootout, "clip4refs"), os.path.join(rootout, "clip5refs"), os.path.join(rootout, "clip6refs"),
outfolders = [os.path.join(rootout, "clip1refs")]

imagelist = []
clip1imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip1.tif"]#, "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"]
imagelist.append(clip1imgs)

searchstrings = ["soy", "corn", "wwheat", "sorghum", "wwheatsoydbl"]
fitmethods = ["SLSQP"]#, "TNC"]

fullpixels = get_px_coords_from_points(clip1imgs[0], "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/clip1_fullcells_points.shp")

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
                    phenological_classificaion(img, outfolder, name + "_" + method, refs, "ENVI", 17, 16, 1000, 0, method, toprint=False, subset=fullpixels)
                else:
                    print "ERROR +++++++++++++++++++++++++++++++++++++++++++++++++++ ERROR"