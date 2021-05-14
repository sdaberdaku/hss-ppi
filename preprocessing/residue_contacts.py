import os
import gzip
import re
from glob import glob
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch
from pathlib import Path

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)

_hydrogen = re.compile("[123 ]*H.*")

def isHydrogen(atm):
    return _hydrogen.match(atm.get_id())


for interaction_distance in [6, 10] :
    for dataset in ["test", "train", "validation"] : # for each dataset partition
        for p_class in ["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"] :
            folder_path = f"data/{dataset}/{p_class}/pdb"
            destination_path = f"data/{dataset}/{p_class}/residue_contacts_{interaction_distance}"
            Path(destination_path).mkdir(parents=True, exist_ok=True)
            for pdb_file in sorted(glob(f"{folder_path}/*.pdb")):
                pdb_id = re.sub("\.pdb$", "", os.path.basename(pdb_file))
                structure = pdb_parser.get_structure(pdb_id, pdb_file)
                print(pdb_id)
                heavy_atoms = [atom for atom in structure.get_atoms() if not isHydrogen(atom)]
                search_tree = NeighborSearch(heavy_atoms)
                residue_nb = {}

                for atom in heavy_atoms:
                    residue = atom.get_parent()
                    if residue not in residue_nb:
                        residue_nb[residue] = set()
                    nb_list = search_tree.search(atom.coord, interaction_distance, "R")
                    residue_nb[residue] = residue_nb[residue].union(nb_list)

                for chain in structure.get_chains():
                    residue_list = [residue for residue in chain.get_residues()]
                    for j in range(len(residue_list)):
                        residue = residue_list[j]
                        if len(residue_list) == 1:
                            seq_residues = [residue]
                        elif 0 < j < len(residue_list)-1:
                            seq_residues = [residue_list[j-1], residue, residue_list[j+1]]
                        elif j == 0:
                            seq_residues = [residue, residue_list[j+1]]
                        else:
                            seq_residues = [residue_list[j-1], residue]
                        residue_nb[residue] = residue_nb[residue].difference(seq_residues)
                
                with gzip.open(f"{destination_path}/{pdb_id}_contacts.txt.gz", "wt") as out_contacts:
                    for residue in structure.get_residues():
                        chain_id = residue.get_full_id()[2]
                        residue_id = residue.get_full_id()[3][1]
                        residue_insertion_code = residue.get_full_id()[3][2]
                        residue_name = residue.get_resname()
                        nb_number = len(residue_nb[residue])
                        out_contacts.write("%s,%s,%s,%s,%d\n" % (chain_id, residue_id, residue_insertion_code, residue_name, nb_number))
