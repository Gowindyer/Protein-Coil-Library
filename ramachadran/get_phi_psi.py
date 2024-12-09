import mdtraj as md
import glob
from tqdm import tqdm
import sys
import numpy as np

_residueName_ = [
    "ALA", "ARG", "ASN", "ASP", "CYS", 
    "GLN", "GLU", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", 
    "SER", "THR", "TRP", "TYR", "VAL"
]
AllResiduePhiPsi = {resn:[] for resn in _residueName_}

PDBFolder = sys.argv[1] 
PDBFileList = glob.glob(PDBFolder + '/*.pdb')
for pdbfile in tqdm(PDBFileList):
    try:
       pdb = md.load(pdbfile)
    except:
        continue
    top = pdb.topology
    atoms = top._atoms
    phi_indices, phi = md.compute_phi(pdb)
    psi_indices, psi = md.compute_psi(pdb)
    idxCalphaPhi = [atoms[phi_idx[2]].index for i, phi_idx in enumerate(phi_indices)]
    idxCalphaPsi = [atoms[psi_idx[1]].index for i, psi_idx in enumerate(psi_indices)]
    resnamePhi = [atoms[phi_idx[2]].residue.name for i, phi_idx in enumerate(phi_indices)]
    resnamePsi = [atoms[psi_idx[1]].residue.name for i, psi_idx in enumerate(psi_indices)]
    coexistIndex = [[idxCalphaPhi.index(idx), idxCalphaPsi.index(idx)] 
                                  for i,idx in enumerate(idxCalphaPhi) 
                                               if idx in idxCalphaPsi
                                               ]
    for idx in coexistIndex:
        resname1 = resnamePhi[idx[0]]
        resname2 = resnamePsi[idx[1]]
        if resname1 in _residueName_:
           AllResiduePhiPsi[resname1].append([phi[0, idx[0]],psi[0,idx[1]]])
np.save('%s/AllResiduePhiPsi.npy'%PDBFolder, AllResiduePhiPsi)
#np.savetxt('%s/AllResiduePhiPsi.dat'%PDBFolder, AllResiduePhiPsi)
