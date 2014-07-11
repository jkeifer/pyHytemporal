import os

############## FUNCTIONS ##############


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











