from utils import find_files
from core import signatureCollection


def get_sigs_in_dir(directory, viname=None, searchstring='mean.ref', recursivesearch=False):
    signatures = signatureCollection(viName=viname)
    sigFiles = find_files(directory, searchstring, recursive=recursivesearch)

    for f in sigFiles:
        signatures.add(f)

    return signatures