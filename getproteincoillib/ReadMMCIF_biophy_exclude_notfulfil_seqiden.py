import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm.app import Modeller
import glob
from tqdm import tqdm
from pdbecif.mmcif_io import CifFileReader
import numpy as np
import json
import os
from Bio import PDB 

def fixPDB(cifentry):
    fixer = PDBFixer(filename=cifentry)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    PDBFile.writeFile(fixer.topology, fixer.positions,open('%s.pdb'%cifentry[:-4], 'w'),keepIds=True)


def getSecondaryStructure(secondaryStructure, sequenceID, chainID,residueSequence,MinFragmentSize):

    # find blank fragment
    blankResidueIndex = np.argwhere(secondaryStructure=='.')
    bsequenceID = sequenceID[blankResidueIndex]   
    bchainID = chainID[blankResidueIndex]
    bresidueSequence = residueSequence[blankResidueIndex]
    
    # find consecutive fragment from sequence ID
    FragmentIndex = np.where(np.diff(bsequenceID[:,0])!= 1)[0]+ 1 
    FsequenceID = np.split(bsequenceID, FragmentIndex) 
    #FSequenceID = [bsequenceID for FragmentIndex]
    FresidueSequence = np.split(bresidueSequence, FragmentIndex)
    FChainID = np.split(bchainID, FragmentIndex)
    

    # get fragment with proper length and delete Pro effect
    residueID = []
    chainID = []
    for i,fresname in enumerate(FresidueSequence):
        if "PRO" in fresname:
           ProIndex = np.argwhere(fresname == 'PRO')[-1][0]
           FsequenceID[i] = FsequenceID[i][ProIndex-1:]
           FChainID[i] = FChainID[i][ProIndex-1:]
           FresidueSequence[i] = FresidueSequence[i][ProIndex-1:]
        if len(FsequenceID[i]) >= MinFragmentSize:
            residueID.append(FsequenceID[i])
            chainID.append(FChainID[i])

    return residueID, chainID

def extract_residues(io,outputPath,residueIndices,chains):
    class SelectByIndex(PDB.Select):
        def accept_residue(self, residue):
            return residue.get_id()[1] in residueIndices  # Use provided indices

        def accept_chain(self, chain):
            return chain.id in chains  # Accept all chains

    io.save(outputPath,SelectByIndex())


# PDB path
folderPDB = sys.argv[1]
# Folder where fragment were saved
cifFilesList = glob.glob(folderPDB + "/*.cif")

outputFolder = 'Fragment_%s'%folderPDB
os.system('mkdir %s'%outputFolder)
resolutionRfree = np.load('%s/structure_resolution_Rfree.npy'%folderPDB, allow_pickle=True).item()

sequenceIdentitymap = np.loadtxt('blastb_fulfilseqidentity/allsequenceIdentity50.dat',dtype=str)
pdbExclution = np.load('blastb_fulfilseqidentity/alldelsequence.npy')
pdbfulfilcondition = sequenceIdentitymap[:,:2].ravel()
faillist = []
for cifentry in tqdm(cifFilesList):
    
    cifCode = cifentry[-8:-4]
    if 'secondary_structure' in cifentry:
        continue
     
    if cifCode not in resolutionRfree:
        continue

    if resolutionRfree[cifCode][1]:
       if resolutionRfree[cifCode][1] > 2:
          continue
    else:
        continue
    if resolutionRfree[cifCode][2]: 
        if resolutionRfree[cifCode][2] > 0.2:
          continue
    elif resolutionRfree[cifCode][3]: 
        if resolutionRfree[cifCode][3] > 0.2:  
          continue
    elif resolutionRfree[cifCode][4]: 
        if resolutionRfree[cifCode][4] > 0.26: 
          continue
    
    # fix PDB
    try:
          fixPDB(cifentry)
    except:
          faillist.append(cifCode + '  Fail to fix PDB\n')
          continue
    
    # generate dssp file in mmcif
    try:
        os.system("mkdssp -v %s/%s.pdb %s_secondary_structure.cif"%(folderPDB, cifCode, cifentry[:-4])) 
    except:
        faillist.append(cifCode + '  Fail to generate dssp\n')
        continue
    
    # read dssp 
    cfr = CifFileReader()
    if not os.path.exists('%s_secondary_structure.cif'%cifentry[:-4]):
        continue
    cif_obj = cfr.read('%s_secondary_structure.cif'%cifentry[:-4], output='cif_wrapper')
    key = list(cif_obj.keys())[0]
    
    try:
       secondaryStructure = np.array(cif_obj[key]._dssp_struct_summary.secondary_structure) 
       sequenceId = np.array(cif_obj[key]._dssp_struct_bridge_pairs.auth_seq_id,dtype=int)
       chainId = np.array(cif_obj[key]._dssp_struct_bridge_pairs.auth_asym_id)
       residueSequence =  np.array(cif_obj[key]._dssp_struct_summary.label_comp_id)
    except:
        faillist.append(cifCode + '  This pdb may be not protein\n')
        continue
    
    if len(sequenceId) == 0:
            faillist.append(cifCode + '  This PDB is not a protein\n')
            continue
    # get blank segment accroding secondary structure
    blankResidueID, blankChainID = getSecondaryStructure(secondaryStructure,
                                                        sequenceId, chainId,
                                                        residueSequence,5) 
    if len(blankResidueID) == 0:
        faillist.append(cifCode + '  There are no proper length of fragments to save\n')
        continue
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', '%s/%s.pdb'%(folderPDB, cifCode))
    io = PDB.PDBIO()
    io.set_structure(structure)
    
    for i,frag in enumerate(blankResidueID):
        pdbchainname = '%s_%s'%(cifCode, blankChainID[i][0][0])
        if pdbchainname not in pdbExclution:
            extract_residues(io,'%s/%s_%02d.pdb'%(outputFolder, cifCode,i), frag, blankChainID[i])

with open('%s/faillist.info'%outputFolder,'w') as fl:
    for line in faillist:
        fl.write(line)

