import pandas as pd
import numpy as np

from glob import glob
import os

folder_path = "./benchmark5/structures/bound"

correlations = []

for receptor_file in sorted(glob("%s/*_r_b.pdb" % folder_path)):
    pdb_id = os.path.basename(receptor_file)[:4]
    print(pdb_id)

    relative_asa = pd.read_csv("./relative_ASA/%s_r_b_rasa.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'ASA'])
    residue_contacts = pd.read_csv("./residue_contacts/%s_r_b_contacts.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'contacts'])

    asa = relative_asa['ASA']
    contacts = residue_contacts['contacts']
    correlations.append(np.corrcoef(asa, contacts)[0, 1])
    print(np.corrcoef(asa, contacts)[0, 1])

for ligand_file in sorted(glob("%s/*_l_b.pdb" % folder_path)):
    pdb_id = os.path.basename(ligand_file)[:4]
    print(pdb_id)

    relative_asa = pd.read_csv("./relative_ASA/%s_l_b_rasa.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'ASA'])
    residue_contacts = pd.read_csv("./residue_contacts/%s_l_b_contacts.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'contacts'])

    asa = relative_asa['ASA']
    contacts = residue_contacts['contacts']
    correlations.append(np.corrcoef(asa, contacts)[0, 1])
    print(np.corrcoef(asa, contacts)[0, 1])


print("\nCorrelation average: %f" % np.mean(correlations))
print("Correlation std: %f" % np.std(correlations))