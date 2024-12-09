from openmm import app
from openmm.app import Modeller
import sys
import numpy as np
from openmm import unit

def gene_psf(psf_file_name,top):
   """
   generate protein structure file for vmd 

   Parameters
   ----------
   psf_file_name: str
      The file name of PSF.

   top: topology
      The topology of OpenMM System.
   """

   atoms = list(top.atoms())
   bonds = list(top.bonds())
   space = ' '
   with open('%s.psf'%psf_file_name,'w') as psf:
         psf.write('PSF\n')
         # ATOM
         #psf.write(space*6)
         psf.write('%8d !NATOM\n'%len(atoms))
         for i in range(len(atoms)):
            psf.write('%8d%7d%7s%s%4s%5s%15d%14d%8d\n'%
                  (atoms[i].index+1,1,atoms[i].residue.name,
                  space*2,atoms[i].name,atoms[i].name,
                  0,100,0))
         psf.write('\n')
         psf.write('%8d !NBOND: bonds\n'%len(bonds))
         for i in range(len(bonds)):
            psf.write('%8d%8d'%(bonds[i][0].index+1,bonds[i][1].index+1))
            if (i+1)%4==0:
               psf.write('\n')
            elif i==int(len(bonds)-1):
               psf.write('\n')

def coarseGrained(pdbFile, pdbname, outputfolder):
    pdb = app.PDBFile(pdbFile)
    topology = pdb.topology
    positions = pdb.getPositions(asNumpy=True)._value
    residues = list(topology.residues())
    modeller = Modeller(pdb.topology, pdb.positions)

    # Get the coordinate of side chain 
    indexRbases = []
    atomRbases = []
    atomH = []
    atomBackbone = []
    residuewcb = []
    massRbases = []
    elementmass = {'N':14.0067, 'C':12.011, 'O':15.999, 'S':32.065 }
    for i_res in residues:
        indrbase = []
        massrbase = []
        for i_a in list(i_res.atoms()):
            if i_a.name not in ['N', 'CA','C', 'O', 'OXT'] and i_a.element.symbol not in ['H', 'D']:
                if i_a.element.symbol in elementmass:
                    i_mass = i_a.element.mass._value
                else:
                    print('%s not in elemantmass for %s.'%(i_a.name,pdbFile))
                indrbase.append(i_a.index)
                massrbase.append(i_mass) 
                residuewcb.append(i_res.index)
                if i_a.name != 'CB':
                   atomRbases.append(i_a)
            if i_a.name in ['N', 'C', 'O', 'OXT']:
                atomBackbone.append(i_a)
            if i_a.element.symbol in ['H', 'D']:
                atomH.append(i_a)
        if len(indrbase) != 0:
           indexRbases.append(indrbase)
           massRbases.append(massrbase)
    if len(indexRbases)!=0:
       coordinateCB = np.zeros((len(indexRbases),3), dtype=float)
       for i, ind in enumerate(indexRbases):
           mass = np.array(massRbases[i]).reshape(-1, 1)
           posi = positions[ind[:],:]
           coordinateCB[i, :] = np.sum(mass*posi, axis=0)/np.sum(mass) 
    allDeleteAtom = atomRbases + atomH + atomBackbone

    # Substitute CB for R-base
    modeller.delete(allDeleteAtom)
    atomCA = [ai for ai in list(modeller.topology.atoms()) if ai.name=='CA']
    for i in range(len(atomCA)-1):
        modeller.topology.addBond(atomCA[i],atomCA[i+1])
    try:
        positions = np.array(modeller.getPositions()._value)
        if len(indexRbases)!= 0:
           indexCB = [ai.index for ai in list(modeller.topology.atoms()) if ai.name == 'CB']
           positions[indexCB,:] = coordinateCB
        app.PDBFile.writeFile(modeller.topology, positions*unit.nanometer, open('%s/%s_cg.pdb'%(outputfolder, pdbname),'w'), keepIds=True)
        gene_psf('%s/%s_cg'%(outputfolder, pdbname), modeller.topology)
    except Exception as e:
        print(pdbname, e)
#pdb = sys.argv[1]
#coarseGrained(pdb, pdb[:-4], '.')
