# computes the sequence adjacency matrix for each pdb file
import os
import re
from glob import glob
from Bio.PDB.PDBParser import PDBParser
import pandas as pd
from pathlib import Path

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)


def is_hetero(res):
    return res.get_full_id()[3][0] != ' '


def residue_name(res):
    return "%s_%d_%s" % (res.get_full_id()[2], res.get_full_id()[3][1], res.get_full_id()[3][2])


for dataset in ["test", "train", "validation"] : # for each dataset partition
    for p_class in ["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"] :
        folder_path = f"data/{dataset}/{p_class}/pdb/"
        destination_path = f"data/{dataset}/{p_class}/sequence_matrices/"
        Path(destination_path).mkdir(parents=True, exist_ok=True)
        for pdb_file in sorted(glob(f"{folder_path}/*.pdb")):
            pdb_id = re.sub("\.pdb$", "", os.path.basename(pdb_file))
            print(pdb_id)
            structure = pdb_parser.get_structure(pdb_id, pdb_file)

            res_names = [residue_name(res) for res in structure.get_residues() if not is_hetero(res)]

            sequence_adjacency_matrix = pd.DataFrame(0, index=res_names, columns=res_names)

            for chain in structure.get_chains():
                residue_list = [res for res in chain.get_residues() if not is_hetero(res)]
                n_residues = len(residue_list)
                for i in range(n_residues-1):
                    res1 = residue_list[i]
                    res2 = residue_list[i+1]
                    n1 = residue_name(res1)
                    n2 = residue_name(res2)
                    sequence_adjacency_matrix.loc[n1, n2] = 1
                    sequence_adjacency_matrix.loc[n2, n1] = 1
            
            print(sequence_adjacency_matrix)
            sequence_adjacency_matrix.to_csv(f"{destination_path}/{pdb_id}_aa_sequence_adjacency_matrix.csv.gz", compression='gzip')
