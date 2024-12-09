from Bio.PDB import MMCIFParser, PPBuilder
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
from tqdm import tqdm
import numpy as np

# Function to extract amino acid sequences from a CIF file
def extract_sequence_from_cif(cif_file):
    parser = MMCIFParser(QUIET=True)  # Quiet mode avoids warnings
    structure = parser.get_structure('protein', cif_file)
    ppb = PPBuilder()  # Polypeptide builder to get sequences
    sequences = []
    
    for chain in structure.get_chains():  # Iterate over chains
        sequence = ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])
        if sequence:
            # Create a SeqRecord for each chain's sequence
            record = SeqRecord(Seq(sequence), id=f"{structure.header['idcode']}_{chain.id}", 
                                description=structure.header['name']+ '; '+structure.header['structure_method'])
            sequences.append(record)
    return sequences

# Directory containing CIF files
cif_directory = sys.argv[1]

# Output FASTA file for extracted sequences
output_fasta = "%s/pdb_sequences_exclude_only_xray_res2_Rp2.fasta"%cif_directory
resolutionRfree = np.load('%s/structure_resolution_Rfree.npy'%cif_directory, allow_pickle=True)
resolutionRfree = resolutionRfree.item()
# Open the output FASTA file
faillist = []
with open(output_fasta, "w") as output_handle:
    # Process each CIF file in the directory
    for cif_file in tqdm(os.listdir(cif_directory)):
        if cif_file.endswith(".cif"):  # Ensure the file is a CIF file
            cifname = cif_file.split('.')[0]
            if cifname not in resolutionRfree:
                continue
            if resolutionRfree[cifname][1]:
               if resolutionRfree[cifname][1] > 2:
                  continue
            else:
                continue
            if resolutionRfree[cifname][2]: 
               if resolutionRfree[cifname][2] > 0.2:
                  continue
            elif resolutionRfree[cifname][3]: 
               if resolutionRfree[cifname][3] > 0.2:  
                  continue
            elif resolutionRfree[cifname][4]: 
               if resolutionRfree[cifname][4] > 0.26: 
                  continue
            cif_path = os.path.join(cif_directory, cif_file)
            try:
               sequences = extract_sequence_from_cif(cif_path)  # Extract sequences
            except Exception as e:
               print(f"Error processing file {cif_file}: {e}")
               faillist.append(cif_file)
               continue
            # Write sequences to FASTA file
            SeqIO.write(sequences, output_handle, "fasta")
with open(cif_directory + '/' + 'failfasta', 'w') as fp:
    fp.write("\n".join(faillist) + "\n")
