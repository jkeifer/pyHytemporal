__author__ = 'phoetrymaster'


import subprocess

nodatain = -3000
nodataout = -3000
inputshape = "'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/DepartmentSelection/ARG_adm/pellegrini.shp'"
inputimg = "'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS 7_2012-2013/argentina_1/test.tif'"
outputimg = "'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS 7_2012-2013/argentina_1/MODIS_pellegrini_clip.tif'"
outprj = "'+proj=utm +zone=20 +datum=WGS84'"
outformat = "ENVI"


#Need to reproject input shapefile to outprj before warping...

subprocess.call("gdalwarp -t_srs {0} -srcnodata {1} -dstnodata {2} -crop_to_cutline -cutline {3} {4} {5} -of {6}".format(outprj, nodatain, nodataout, inputshape, inputimg, outputimg, outformat), shell=True)