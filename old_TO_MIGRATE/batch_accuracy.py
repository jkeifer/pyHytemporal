import os
import sys
from datetime import datetime as dt
import accuracy_assessment
import multiprocessing

def main():

    start = dt.now()
    print(start)

    maxworkers = 5

    nodata = -3000
    searchstringsvals = [("soy.", 5), ("wwheat.", 24), ("corn", 1)]
    outputlocation = r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/ALLCLASSIFIEDIMAGES"

    searchdirs = [r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip1refs/KansasNDVI_2012_clip6_SLSQP'
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip2refs/KansasNDVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip3refs/KansasNDVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip4refs/KansasNDVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip5refs/KansasNDVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/clip6refs/KansasNDVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasEVI_2012_clip6_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip1_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip2_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip3_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip4_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip5_SLSQP',
                  r'/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/MODIS_KANSAS_2007-2012/reprojected/Classified/meanrefs/KansasNDVI_2012_clip6_SLSQP']

    cropimages = [r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip1.tif",
                  r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip2.tif",
                  r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip3.tif",
                  r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip4.tif",
                  r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip5.tif",
                  r"/Users/phoetrymaster/Documents/School/Geography/Thesis/Data/polygonclip_20130929223024_325071991/resampled/newclips/2012clip6.tif"]

    searchdirs.reverse()

    while len(searchdirs) > 0:
        processes = []
        for i in range(0, maxworkers):
            searchdir = searchdirs.pop()

            if "_clip1_" in searchdir:
                cropimg = cropimages[0]
            if "_clip2_" in searchdir:
                cropimg = cropimages[1]
            if "_clip3_" in searchdir:
                cropimg = cropimages[2]
            if "_clip4_" in searchdir:
                cropimg = cropimages[3]
            if "_clip5_" in searchdir:
                cropimg = cropimages[4]
            if "_clip6_" in searchdir:
                cropimg = cropimages[5]

            p = multiprocessing.Process(target=accuracy_assessment.main, args=(searchdir, cropimg, searchstringsvals, nodata))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    end = dt.now() - start
    print(end)


if __name__ == '__main__':
    sys.exit(main())
