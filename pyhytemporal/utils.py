import os

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


def create_output_dir(root, name):
    dirpath = os.path.join(root, name)

    if os.path.isdir(dirpath):
        count = 1
        dirpath_ = dirpath + "_"
        while 1:
            dirpath = dirpath_ + str(count)
            count += 1
            if not os.path.isdir(dirpath):
                break

    os.makedirs(dirpath)
    return dirpath