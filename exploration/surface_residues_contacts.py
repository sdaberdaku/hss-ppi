import pandas as pd
import numpy as np

from glob import glob
import os
import sys

folder_path = "./benchmark5/structures/bound"

out_path = "./"

for receptor_file in sorted(glob("%s/*_r_b.pdb" % folder_path)):
    pdb_id = os.path.basename(receptor_file)[:4]
    print(pdb_id)

    relative_asa = pd.read_csv("./relative_ASA/%s_r_b_rasa.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'ASA'])
    residue_contacts = pd.read_csv("./residue_contacts/%s_r_b_contacts.txt" % pdb_id, sep=",", header=0, names=['chain', 'seq', 'ins', 'resname', 'contacts'])

    asa = relative_asa['ASA']
    contacts = residue_contacts['contacts']

    mean_relative_asa = np.mean(asa)
    std_relative_asa = np.std(asa)

    mean_contacts = np.mean(contacts)
    std_contacts = np.std(contacts)

    surface_res_high_contacts = []

    for i in range(relative_asa.shape[0]):
        # print(relative_asa['ASA'][i])
        # print(residue_contacts['contacts'][i])
        if (relative_asa['ASA'][i] < mean_relative_asa) and (residue_contacts['contacts'][i] >= mean_contacts):
            surface_res_high_contacts.append(i)

    print(surface_res_high_contacts)

