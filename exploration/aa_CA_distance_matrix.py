# computes the truth value for every residue in the pdb files
import os
import re
from glob import glob
from Bio.PDB.PDBParser import PDBParser
import pandas as pd
import math

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)

_hydrogen = re.compile("[123 ]*H.*")


def isHydrogen(atm):
    return _hydrogen.match(atm.get_id())


def isHetero(res):
    return res.get_full_id()[3][0] != ' '

def residue_name(res):
    return "%s_%d_%s" % (res.get_full_id()[2], res.get_full_id()[3][1], res.get_full_id()[3][2])

folder_path = "./benchmark5/structures/bound"
destination_path = "./ca_distance_matrices"

for pdb_file in sorted(glob("%s/*_b.pdb" % folder_path)):
    pdb_id = os.path.basename(pdb_file)[:6]
    print(pdb_id)
    structure = pdb_parser.get_structure(pdb_id, pdb_file)

    res_names = [residue_name(res) for res in structure.get_residues() if not isHetero(res)]

    distance_matrix = pd.DataFrame(0, index=res_names, columns=res_names)

    residue_list = [res for res in structure.get_residues() if not isHetero(res)]

    n_residues = len(residue_list)

    for i in range(n_residues-1):
        for j in range(i+1, n_residues):
            res1 = residue_list[i]
            res2 = residue_list[j]
            n1 = residue_name(res1)
            n2 = residue_name(res2)

            min_distance = math.inf
            try:
                distance = res1['CA'] - res2['CA']
                if distance < min_distance:
                    min_distance = distance
            except KeyError:
                # no CA atom in one of the two
                if 'CA' in res1:
                    heavy_atoms_1 = [res1['CA']]
                else:
                    print("missing CA in residue:")
                    print(res1.get_full_id())
                    heavy_atoms_1 = [atom for atom in res1.get_atoms() if not isHydrogen(atom)]
                if 'CA' in res2:
                    heavy_atoms_2 = [res2['CA']]
                else:
                    print("missing CA in residue:")
                    print(res2.get_full_id())
                    heavy_atoms_2 = [atom for atom in res2.get_atoms() if not isHydrogen(atom)]
                for atm1 in heavy_atoms_1:
                    for atm2 in heavy_atoms_2:
                        distance = atm1 - atm2
                        if distance < min_distance:
                            min_distance = distance

            distance_matrix.loc[n1, n2] = min_distance
            distance_matrix.loc[n2, n1] = min_distance

    print(distance_matrix)
    distance_matrix.to_csv("%s/%s_aa_distance_matrix.csv" % (destination_path, pdb_id))

