import numpy as np
from glob import glob 
from tqdm import tqdm
import sys
from coarseGrainedPDB import coarseGrained
import os

folderPDB = sys.argv[1]
pdbFileList = glob(folderPDB + "/*.pdb")
outputFolder = '%s_CG'%folderPDB

os.system('mkdir %s'%outputFolder)
for pdb in tqdm(pdbFileList):
    pdbname = pdb.split('/')[2]
    pdbname = pdbname.split('.')[0]
    coarseGrained(pdb, pdbname, outputFolder)
#coarseGrained('6M9Z_01.pdb', '6M9Z_01', '.')
