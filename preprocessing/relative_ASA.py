import re
import gzip
import os
from glob import glob
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from pathlib import Path
from itertools import zip_longest

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)


def is_hetero(res):
    return res.get_full_id()[3][0] != ' '


for dataset in ["test", "train", "validation"] : # for each dataset partition
    for p_class in ["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"] :
        folder_path = f"data/{dataset}/{p_class}/pdb"
        destination_path = f"data/{dataset}/{p_class}/relative_ASA"
        Path(destination_path).mkdir(parents=True, exist_ok=True)
        for pdb_file in sorted(glob(f"{folder_path}/*.pdb")):
            pdb_id = re.sub("\.pdb$", "", os.path.basename(pdb_file))
            structure = pdb_parser.get_structure(pdb_id, pdb_file)
            dssp = DSSP(structure[0], pdb_file, acc_array = "Wilke")

            with gzip.open(f"{destination_path}/{pdb_id}_rasa.txt.gz", "wt") as out_rasa:
                for dssp_key, residue in zip_longest(dssp.keys(), structure.get_residues()):
                    chain_id = residue.get_full_id()[2]
                    residue_id = residue.get_full_id()[3][1]
                    residue_insertion_code = residue.get_full_id()[3][2]
                    residue_name = residue.get_resname()
                    if dssp_key is None:
                        relative_ASA = -1
                    else:
                        dssp_data = dssp[dssp_key]
                        relative_ASA = dssp_data[3]
                        if relative_ASA == 'NA':
                            relative_ASA = -1
                    out_rasa.write("%s,%s,%s,%s,%s\n" % (chain_id, residue_id, residue_insertion_code, residue_name, relative_ASA))
