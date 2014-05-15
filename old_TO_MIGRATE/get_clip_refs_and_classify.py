from get_ref_values_v2 import *
from create_rule_image_multiprocessed_bypx import read_reference_file, get_sort_dates_values, phenological_classificaion
import os
from pyhytemporal.utils import find_files

outroot = "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Refs/2012"

#clip1
wwheat = [(39, 37), (45, 50), (11, 54), (4, 19), (75, 34)]
corn = [(48, 44), (10, 68), (8, 37), (67, 66), (83, 21)]
soy = [(42, 45), (31, 86), (7, 64), (55, 83), (73, 73)]
wwheatsoydbl = [(32, 62), (38, 55), (25, 56), (39, 30), (97, 3)]
sorghum = [(82, 52), (70, 38), (49, 33), (36, 28), (28, 27)]

clip1locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip1imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip1.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip1.tif"]
clip1out = os.path.join(outroot, "clip1")


#clip2
wwheat = [(2, 11), (6, 28), (63, 29), (23, 55), (6, 55), (3, 36)]
corn = [(3, 36), (79, 53), (50, 50), (66, 25), (70, 4), (65, 63)]
soy = [(6, 25), (6, 14), (97, 35), (85, 49), (60, 25)]
wwheatsoydbl = [(41, 53), (73, 61), (44, 32), (47, 64)]
sorghum = []

clip2locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip2imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip2.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip2.tif"]
clip2out = os.path.join(outroot, "clip2")


#clip3
wwheat = [(9, 28), (12, 22), (21, 61), (15, 43)]
corn = [(52, 89), (30, 65), (16, 67), (57, 40), (11, 41)]
soy =[(30, 68), (14, 72), (18, 40), (18, 36), (88, 43)]
wwheatsoydbl = [(61, 82), (50, 57), (67, 39), (37, 62), (43, 58)]
sorghum = []

clip3locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip3imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip3.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip3.tif"]
clip3out = os.path.join(outroot, "clip3")


#clip4
wwheat = [(46, 63), (42, 74), (26, 74), (67, 57), (69, 86), (76, 42)]
corn = [(60, 71), (56, 71), (53, 71), (46, 68), (46, 57)]
soy = [(56, 82), (52, 82), (35, 78), (32, 74), (18, 70)]
wwheatsoydbl = [(39, 70), (39, 67), (26, 53), (52, 92), (2, 55)]
sorghum = [(56, 68), (5, 62), (10, 42), (6, 51), (10, 18)]

clip4locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip4imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip4.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip4.tif"]
clip4out = os.path.join(outroot, "clip4")


#clip5
wwheat = [(52, 82), (51, 40), (22, 67), (14, 76), (8, 62), (24, 91)]
corn = [(35, 67), (32, 57), (65, 44), (20, 35), (86, 81)]
soy = [(24, 32), (23, 35), (20, 39), (48, 23), (91, 59)]
wwheatsoydbl = []
sorghum = [(38, 74), (43, 61), (99, 78), (3, 80), (44, 95), (56, 76)]

clip5locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip5imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip5.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip5.tif"]
clip5out = os.path.join(outroot, "clip5")


#clip6
wwheat = [(59, 78), (18, 66), (39, 28), (3, 17), (36, 11)]
corn = [(70, 70), (38, 62), (45, 32), (46, 21), (60, 18)]
soy = [(16, 60), (39, 55), (70, 85), (48, 28), (56, 14)]
wwheatsoydbl = [(31, 33), (7, 87), (4, 94), (90, 11)]
sorghum = [(77, 43), (10, 94), (32, 70), (33, 18), (40, 10)]

clip6locs = {"wwheat": wwheat, "corn": corn, "soy": soy, "wwheatsoydbl": wwheatsoydbl, "sorghum": sorghum}
clip6imgs = ["/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasNDVI_2012_clip6.tif", "/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/clips/KansasEVI_2012_clip6.tif"]
clip6out = os.path.join(outroot, "clip6")


clips = [(clip1locs, clip1imgs, clip1out), (clip2locs, clip2imgs, clip2out), (clip3locs, clip3imgs, clip3out), (clip4locs, clip4imgs, clip4out), (clip5locs, clip5imgs, clip5out), (clip6locs, clip6imgs, clip6out)]

for locs, imgs, out in clips:
    for img in imgs:
        print "Processing {0}...".format(img)
        if "NDVI" in img:
            postfix = "_NDVI"
        else:
            postfix = "_EVI"
        get_reference_curves(img, locs, 17, 16, outdir=out, filepostfix=postfix)

means = find_files(outroot, "mean.ref")
searchstrings = (["soy", "wwheat", "sorghum", "corn", "wwheatsoydbl"], ["_EVI", "_NDVI"])

for string in searchstrings[0]:
    for vi in searchstrings[1]:
        cropmeans = []
        for mean in means:
            if string + vi in mean:
                cropmeans.append(mean)
        cropname = string + vi
        refvals = []
        for cropmean in cropmeans:
            error, vals = read_reference_file(cropmean)
            if not error:
                refvals.append(vals)
        listedvals = []
        for ref in refvals:
            x, y = get_sort_dates_values(ref)
            listedvals.append(y)
        print cropname
        write_mean_ref_to_txt(cropname, listedvals, 17, 16, outroot, comment="Generated by clipref.py version 0.1")
