import sys
from pdbecif.mmcif_io import CifFileReader
import numpy as np
import glob
from tqdm import tqdm

def strtofloat(element):
    try:
        element = float(element[0])
    except:
        element = None
    return element

cifFolder = sys.argv[1]
cifFilesList = glob.glob(cifFolder + "/*.cif")

resolution_and_Rfree = {}
for ciffile in tqdm(cifFilesList):

    cifRead = CifFileReader()
    cifObj = cifRead.read(ciffile, output='cif_wrapper')
    key = list(cifObj.keys())[0]
    exptlmethod = cifObj[key]._exptl.method
    strexptlmethod = ''.join(exptlmethod)
    if 'DIFFRACTION' in strexptlmethod:
        try:
            technics = cifObj[key]._refine.pdbx_refine_id
            resolution = strtofloat(cifObj[key]._refine.ls_d_res_high)
            Rfactorobs = strtofloat(cifObj[key]._refine.ls_R_factor_obs)
            Rfactorwork = strtofloat(cifObj[key]._refine.ls_R_factor_R_work)
            Rfactorfree = strtofloat(cifObj[key]._refine.ls_R_factor_R_free)
            entryID = cifObj[key]._refine.entry_id
            resolution_and_Rfree[entryID[0]] = [technics[0], resolution, Rfactorobs, 
                                                Rfactorwork, Rfactorfree] 
        except:
            print(ciffile)
np.save('%s/structure_resolution_Rfree.npy'%cifFolder, resolution_and_Rfree)
