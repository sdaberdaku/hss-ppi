import pandas as pd
from glob import glob
import os
import re
import matplotlib.pyplot as plt
from utils.load_protein_data import graph_weights
import time
import numpy as np
from scipy.sparse.linalg.eigen.arpack import eigsh

# https://en.omegaxyz.com/2021/01/25/von-Neumann-Graph-Entropy-Python-Implementation/
def normalized_laplacian(adj_matrix):
    nodes_degree = np.sum(adj_matrix, axis=1)
    nodes_degree_sqrt = 1/np.sqrt(nodes_degree)
    degree_matrix = np.diag(nodes_degree_sqrt)
    eye_matrix = np.eye(adj_matrix.shape[0])
    return eye_matrix - degree_matrix * adj_matrix * degree_matrix


def unnormalized_laplacian(adj_matrix):
    nodes_degree = np.sum(adj_matrix, axis=1)
    degree_matrix = np.diag(nodes_degree)
    return degree_matrix - adj_matrix


def VNGE_exact(adj_matrix):
    start = time.time()
    nodes_degree = np.sum(adj_matrix, axis=1)
    c = 1.0 / np.sum(nodes_degree)
    laplacian_matrix = c * unnormalized_laplacian(adj_matrix)
    eigenvalues, _ = np.linalg.eig(laplacian_matrix)
    eigenvalues[eigenvalues < 0] = 0
    pos = eigenvalues > 0
    H_vn = - np.sum(eigenvalues[pos] * np.log2(eigenvalues[pos]))
    print('H_vn exact:', H_vn)
    print('Time:', time.time() - start)
    return H_vn


def VNGE_FINGER(adj_matrix):
    start = time.time()
    nodes_degree = np.sum(adj_matrix, axis=1)
    c = 1.0 / np.sum(nodes_degree)
    edge_weights = 1.0 * adj_matrix[np.nonzero(adj_matrix)]
    approx = 1.0 - np.square(c) * (np.sum(np.square(nodes_degree)) + np.sum(np.square(edge_weights)))
    laplacian_matrix = unnormalized_laplacian(adj_matrix)
    '''
    eigenvalues, _ = np.linalg.eig(laplacian_matrix)  # the biggest reduction
    eig_max = c * max(eigenvalues)
    '''
    eig_max, _ = eigsh(laplacian_matrix, 1, which='LM')
    eig_max = eig_max[0] * c
    H_vn = - approx * np.log2(eig_max)
    print('H_vn approx:', H_vn)
    print('Time:', time.time() - start)





datasets = ["train", "validation", "test"]
protein_class = "ER"
is_receptor = True
is_bound = True
representation_type = "hierarchical"
distance_threshold = 6
lr = "r" if is_receptor else "l"
bu = "b" if is_bound else "u"

entropies = []
for dataset in datasets :
    for pdb_file in sorted(glob(f"data/{dataset}/{protein_class}/pdb/*_{lr}_{bu}.pdb")):
        pdb_id = re.sub("\.pdb$", "", os.path.basename(pdb_file))
        print(f"Loading {pdb_id}...")
        dist_matrix = pd.read_csv(f"data/{dataset}/{protein_class}/distance_matrices/{pdb_id}_aa_distance_matrix.csv.gz",
                                  header=0, index_col=0, compression='gzip')
        seq_matrix = pd.read_csv(f"data/{dataset}/{protein_class}/sequence_matrices/{pdb_id}_aa_sequence_adjacency_matrix.csv.gz",
                                 header=0, index_col=0, compression='gzip')

        adj_m = graph_weights(dist_matrix, seq_matrix, representation_type, distance_threshold)

        print(adj_m.shape)

        e = VNGE_exact(adj_m)
        entropies.append(e)

print(f"mean: {np.mean(entropies)}")
print(f"variance: {np.var(entropies)}")
