# computes the aa distance matrix for every residue pair in each pdb file
import os
import re
from glob import glob
from Bio.PDB.PDBParser import PDBParser
import pandas as pd
import math
from pathlib import Path
from multiprocessing import Pool, cpu_count

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)

_hydrogen = re.compile("[123 ]*H.*")


def is_hydrogen(atm):
    return _hydrogen.match(atm.get_id())


def is_hetero(res):
    return res.get_full_id()[3][0] != ' '


def residue_name(res):
    return "%s_%d_%s" % (res.get_full_id()[2], res.get_full_id()[3][1], res.get_full_id()[3][2])


def compute_distance_matrix(pdb_file, pdb_id, dest_path):
    # print(pdb_id)
    structure = pdb_parser.get_structure(pdb_id, pdb_file)
    # we filter out HETATM entries since they are not standard residue atoms
    residue_list = [res for res in structure.get_residues() if not is_hetero(res)]

    res_names = [residue_name(res) for res in residue_list]

    distance_matrix = pd.DataFrame(0, index=res_names, columns=res_names)

    n_residues = len(residue_list)
    # precompute the list of heavy atoms (non-hydrogens) for each AA
    residue_atoms = {}
    for i in range(n_residues) :
        res = residue_list[i]
        heavy_atoms = [atom for atom in res.get_atoms() if not is_hydrogen(atom)]
        residue_atoms[i] = heavy_atoms

    for i in range(n_residues - 1):
        res1 = residue_list[i]
        n1 = residue_name(res1)
        heavy_atoms_1 = residue_atoms[i]
        for j in range(i + 1, n_residues):
            res2 = residue_list[j]
            n2 = residue_name(res2)
            heavy_atoms_2 = residue_atoms[j]
            min_distance = math.inf
            for atm1 in heavy_atoms_1:
                for atm2 in heavy_atoms_2:
                    distance = atm1 - atm2
                    if distance < min_distance:
                        min_distance = distance

            distance_matrix.loc[n1, n2] = min_distance
            distance_matrix.loc[n2, n1] = min_distance

    # print(distance_matrix)
    distance_matrix.to_csv(f"{dest_path}/{pdb_id}_aa_distance_matrix.csv.gz", compression='gzip')


def wrapped_function_call(args):
    """
    we need to wrap the call to unpack the parameters
    we build before as a tuple for being able to use pool.map
    """
    compute_distance_matrix(*args)


pool = Pool(cpu_count())  # machine core number

all_args = []
for dataset in ["test", "train", "validation"] : # for each dataset partition
    for p_class in ["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"] :
        folder_path = f"data/{dataset}/{p_class}/pdb/"
        destination_path = f"data/{dataset}/{p_class}/distance_matrices/"
        Path(destination_path).mkdir(parents=True, exist_ok=True)
        all_args += [(pdb_file, re.sub("\.pdb$", "", os.path.basename(pdb_file)), destination_path) for pdb_file in glob(f"{folder_path}/*.pdb")]

print(all_args)
pool.map(wrapped_function_call, all_args)

