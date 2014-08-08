import os
from constants import *

class log(object):
    """
    Logging utility. Creates a .txt file at init at supplied path. Default verbosity is ERROR (3) which will only allow
    messages of same or higher level display on screen. All messages, regardless of level, will be written to the .txt
    file.
    """

    def __init__(self, txtpath, verbosity=ERROR):
        self.file = open(txtpath, 'w')
        self.verbosity = verbosity

    def __del__(self):
        self.file.close()

    def log(self, message, level=INFO, nolinebreak=False):
        if nolinebreak:
            linebreak = ""
        else:
            linebreak = "\n"
        if level >= self.verbosity:
            print(message)
        self.file.write(message + linebreak)


def find_files(searchdir, ext, recursive=True):
    foundfiles = []

    if recursive:
        for root, dirs, files in os.walk(searchdir):
            for f in files:
                if f.upper().endswith(ext.upper()):
                    foundfile = os.path.join(root, f)
                    foundfiles.append(foundfile)
    else:
        files = os.listdir(searchdir)
        for f in files:
            if f.upper().endswith(ext.upper()):
                foundfile = os.path.join(searchdir, f)
                foundfiles.append(foundfile)

    return foundfiles


def unique_name(root, name, ext="", usetime=False):

    if usetime:
        import datetime
        name = name + datetime.datetime.now().strftime("_%Y-%m-%d_%H%M")

    fullpath = os.path.join(root, name + ext)

    if os.path.exists(fullpath):
        count = 1
        name_ = name + "_"
        while 1:
            fullpath = os.path.join(root, name_ + str(count) + ext)
            count += 1
            if not os.path.exists(fullpath):
                break
    else:
        pass

    return fullpath


def create_output_dir(root, name, usetime=False):
    dirpath = unique_name(root, name, usetime=usetime)
    os.makedirs(dirpath)
    return dirpath


def band_number_to_doy(bandnumber, startDOY, imageryinterval):
    """
    """
    calcdoy = (bandnumber - 1) * imageryinterval + startDOY

    if calcdoy > 0:
        calcdoy -= calcdoy % 365 % imageryinterval - 1

    return calcdoy


def change_geotransform(originaltransform, newxmin, newymin):
    """

    """
    #TODO Docstring

    newtransform = list(originaltransform)

    newtransform[0] += newtransform[1] * newxmin
    newtransform[3] += newtransform[5] * newymin

    return tuple(newtransform)