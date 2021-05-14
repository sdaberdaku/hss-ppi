# computes the truth value for every residue in the pdb files
import os
import gzip
import re
from glob import glob
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch
from pathlib import Path

pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)
# regex used to distinguish hydrogen atoms
_hydrogen = re.compile("[123 ]*H.*")

def is_hydrogen(atm):
    return _hydrogen.match(atm.get_id())

def is_hetero(res):
    return res.get_full_id()[3][0] != ' '

def get_binding_partner_id(pdb_id) :
    partner_id = pdb_id.split("_")
    partner_id[1] = "l" if partner_id[1] is "r" else "r"
    return "_".join(partner_id)

interaction_distance = 6.0

for dataset in ["test", "train", "validation"] : # for each dataset partition
    for p_class in ["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"] :
        folder_path = f"data/{dataset}/{p_class}/pdb/"
        for pdb_file_1 in sorted(glob(f"{folder_path}/*.pdb")):
            pdb_id_1 = re.sub("\.pdb$", "", os.path.basename(pdb_file_1))
            pdb_id_2 = get_binding_partner_id(pdb_id_1)
            print(f"Computing interacting interface for {pdb_id_1} - {pdb_id_2}.")
            pdb_file_2 = f"data/benchmark5/structures/{pdb_id_2}.pdb"

            p1_structure = pdb_parser.get_structure(pdb_id_1, pdb_file_1)
            p2_structure = pdb_parser.get_structure(pdb_id_2, pdb_file_2)

            p1_heavy_atoms = [atom for atom in p1_structure.get_atoms() if not is_hydrogen(atom)]
            p2_heavy_atoms = [atom for atom in p2_structure.get_atoms() if not is_hydrogen(atom)]

            p1_tree = NeighborSearch(p1_heavy_atoms)
            p1_interface = set()

            for atom in p2_heavy_atoms:
                nb_list = p1_tree.search(atom.coord, interaction_distance, "R")
                p1_interface = p1_interface.union(nb_list)

            if not p1_interface:
                print(f"No interface residues for PDB id {pdb_id} in {p_class} class, {dataset} dataset!")

            destination_path = f"data/{dataset}/{p_class}/interaction_interface"
            Path(destination_path).mkdir(parents=True, exist_ok=True)
            with gzip.open(f"{destination_path}/{pdb_id_1}_interface.txt.gz", "wt") as out_p1:
                for residue in p1_structure.get_residues():
                    if is_hetero(residue):
                        continue
                    chain_id = residue.get_full_id()[2]
                    residue_id = residue.get_full_id()[3][1]
                    residue_insertion_code = residue.get_full_id()[3][2]
                    residue_name = residue.get_resname()
                    is_interface = 1 if residue in p1_interface else 0
                    out_p1.write("%s,%s,%s,%s,%d\n" % (chain_id, residue_id, residue_insertion_code, residue_name, is_interface))

