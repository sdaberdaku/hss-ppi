import scipy.sparse as sp
import numpy as np
import pandas as pd
from glob import glob
import os
import re
import math
from sklearn.preprocessing import StandardScaler


def graph_weights(distance_matrix, sequence_matrix, representation_type, threshold):
    if (representation_type == "sequence") :
        weights_matrix = distance_matrix.applymap(lambda x: 0)
    elif (representation_type == "contact_map") : 
        weights_matrix = distance_matrix.applymap(lambda x: 1 if  x <= threshold else 0)
    else:
        weights_matrix = distance_matrix.applymap(lambda x: 1/(1+x) if  x <= threshold else 0)
    weights_matrix = weights_matrix.where(sequence_matrix == 0, 1)
    return weights_matrix


def sample_mask(idx, l):
    """
    Create mask.
    """
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)


def load_protein_data(aa_index_file, representation_type, protein_class, is_bound, is_receptor, distance_threshold):
    # We first read the AAIndex 1 Features from file.
    aaindex_features = pd.read_csv(aa_index_file, header=0, index_col=0)
    # ASX is a notation used when it is unclear whether the current amino acid is ASN or ASP. 
    # GLX is a notation used when it is unclear whether the current amino acid is GLN or GLU. 
    # In these cases, we simply assign the mean of the corresponding feature values.
    aaindex_features["ASX"] = aaindex_features.loc[:, "ASN":"ASP"].mean(axis=1)
    aaindex_features["GLX"] = aaindex_features.loc[:, "GLN":"GLU"].mean(axis=1)


    with open("data/selected_features_2.txt", "r") as in_features:
        selected_features = [x.strip() for x in in_features.readlines()]
    # We are only interested in a subset of all the features of the AAIndex.
    aaindex_features = aaindex_features.loc[selected_features, :]
    # We also add the one-hot-encoding to the selected features
    one_hot_encoding = pd.read_csv("data/one_hot_encoding.csv", header=0, index_col=0)
    features_df = aaindex_features.append(one_hot_encoding)
    
    # feature_names = features_df.index
    print(aaindex_features)

    labels = np.empty(shape=(0, 1))
    # Since we are creating one big block-diagonal matrix representing all the proteins in the dataset,
    # we need a way to identify to which protein does each amino acid belong to. We use the groups 
    # variable for this scope, for instance groups["Training_Set"][7] will contain the name of the 
    # protein of the training set which contains the amino acid at index 7.  
    groups = {"train": [],
              "validation": [],
              "test": []}
    adj_matrix = {"train": None,
                  "validation": None,
                  "test": None}
    features = {"train": None,
                "validation": None,
                "test": None}
    lr = "r" if is_receptor else "l"
    bu = "b" if is_bound else "u"
    # Here we import all the proteins from the training, validation and test sets.
    for dataset in ["train", "validation", "test"]:
        adj_list = []
        features_list = []
        for pdb_file in sorted(glob(f"data/{dataset}/{protein_class}/pdb/*_{lr}_{bu}.pdb")):
            pdb_id = re.sub("\.pdb$", "", os.path.basename(pdb_file))
            print(f"Loading {pdb_id}...")
            seq_matrix = pd.read_csv(f"data/{dataset}/{protein_class}/sequence_matrices/{pdb_id}_aa_sequence_adjacency_matrix.csv.gz", header=0, index_col=0, compression='gzip')
            dist_matrix = pd.read_csv(f"data/{dataset}/{protein_class}/distance_matrices/{pdb_id}_aa_distance_matrix.csv.gz", header=0, index_col=0, compression='gzip')
            interface = pd.read_csv(f"data/{dataset}/{protein_class}/interaction_interface/{pdb_id}_interface.txt.gz", names=["chain", "seq", "ins", "name", "interface"], compression='gzip')
            relative_ASA = pd.read_csv(f"data/{dataset}/{protein_class}/relative_ASA/{pdb_id}_rasa.txt.gz", names=["chain", "seq", "ins", "name", "rASA"], compression='gzip')
            ASA = pd.read_csv(f"data/{dataset}/{protein_class}/ASA/{pdb_id}_asa.txt.gz", names=["chain", "seq", "ins", "name", "ASA"], compression='gzip')
            PHI = pd.read_csv(f"data/{dataset}/{protein_class}/PHI/{pdb_id}_phi.txt.gz", names=["chain", "seq", "ins", "name", "phi"], compression='gzip')
            PSI = pd.read_csv(f"data/{dataset}/{protein_class}/PSI/{pdb_id}_psi.txt.gz", names=["chain", "seq", "ins", "name", "psi"], compression='gzip')
            if (distance_threshold >= 10) :
                residue_contacts = pd.read_csv(f"data/{dataset}/{protein_class}/residue_contacts_10/{pdb_id}_contacts.txt.gz", names=["chain", "seq", "ins", "name", "contacts"], compression='gzip')
            else : 
                residue_contacts = pd.read_csv(f"data/{dataset}/{protein_class}/residue_contacts_6/{pdb_id}_contacts.txt.gz", names=["chain", "seq", "ins", "name", "contacts"], compression='gzip')
            current_features_list = []
    
            for _, row in interface.iterrows():
                current_features_list.append(features_df.loc[:, row["name"]])
            
            current_features = pd.concat(current_features_list, ignore_index=True, axis=1).T
            
            current_features["psi"] = PSI["psi"]
            current_features["phi"] = PHI["phi"]
            current_features["contacts"] = residue_contacts["contacts"]
            current_features["ASA"] = ASA["ASA"]
            current_features["rASA"] = relative_ASA["rASA"]

            current_features = sp.csr_matrix(current_features)
            features_list.append(current_features)
            
            n_nodes = interface.shape[0]
            groups[dataset] += [pdb_id] * n_nodes

            labels = np.append(labels, [[int(x)] for x in interface.loc[:, "interface"].values], axis=0)

            current_adj_matrix = graph_weights(dist_matrix, seq_matrix, representation_type, distance_threshold)
            current_adj_matrix = sp.csr_matrix(current_adj_matrix.values)
            adj_list.append(current_adj_matrix)

        features[dataset] = sp.vstack(features_list)
        adj_matrix[dataset] = sp.block_diag(adj_list, format="csr")

    n_samples_train = features["train"].shape[0]
    n_samples_validation = features["validation"].shape[0]
    n_samples_test = features["test"].shape[0]

    adj_matrix = sp.block_diag([adj_matrix["train"],
                                adj_matrix["validation"],
                                adj_matrix["test"]], format="csr")

    # Standardize the training, validation and test sets
    std_scaler = StandardScaler()
    features["train"] = std_scaler.fit_transform(features["train"].todense())
    features["validation"] = std_scaler.transform(features["validation"].todense())
    features["test"] = std_scaler.transform(features["test"].todense())

    features = sp.vstack([sp.csc_matrix(features["train"]),
                          sp.csc_matrix(features["validation"]),
                          sp.csc_matrix(features["test"])]).tolil()

    # Here we determine the boundaries of each set.
    idx_train = range(n_samples_train)
    idx_val = range(n_samples_train, n_samples_train + n_samples_validation)
    idx_test = range(n_samples_train + n_samples_validation, n_samples_train + n_samples_validation + n_samples_test)
    # Then we build bit-masks for each set
    train_mask = sample_mask(idx_train, labels.shape[0])
    val_mask = sample_mask(idx_val, labels.shape[0])
    test_mask = sample_mask(idx_test, labels.shape[0])
    # Finally, we prepare the labels for each set.
    y_train = np.zeros(labels.shape)
    y_val = np.zeros(labels.shape)
    y_test = np.zeros(labels.shape)
    y_train[train_mask, :] = labels[train_mask, :]
    y_val[val_mask, :] = labels[val_mask, :]
    y_test[test_mask, :] = labels[test_mask, :]

    return adj_matrix, features, y_train, y_val, y_test, train_mask, val_mask, test_mask, groups["train"], groups["validation"], groups["test"]
