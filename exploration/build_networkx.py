import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import os


def graph_weights(distance_matrix, sequence_matrix, threshold):
    weights_matrix = distance_matrix.applymap(lambda x: 1 / (1 + x) if x <= threshold else 0)
    weights_matrix = weights_matrix.where(sequence_matrix == 0, 1)
    return weights_matrix


with open("selected_features.txt", "r") as in_features:
    feature_list = in_features.read().splitlines()

print(feature_list)

aaindex_features = pd.read_csv("AAIndexFeatures.csv", header=0, index_col=0)

aaindex_features = aaindex_features.loc[feature_list, :]

# Sometimes it is not possible two differentiate two closely related amino acids, therefore we have the special cases:
# asparagine/aspartic acid - ASX - B
# glutamine/glutamic acid - GLX - Z
aaindex_features["ASX"] = aaindex_features.loc[:, "ASN":"ASP"].mean(axis=1)
aaindex_features["GLX"] = aaindex_features.loc[:, "GLN":"GLU"].mean(axis=1)

print(aaindex_features)

folder_path = "./benchmark5/structures/bound"
destination_path = "./networkx"

for pdb_file in sorted(glob("%s/*_b.pdb" % folder_path)):
    pdb_id = os.path.basename(pdb_file)[:6]
    print(pdb_id)

    seq_matrix = pd.read_csv("./sequence_matrices/%s_aa_sequence_adjacency_matrix.csv" % pdb_id, header=0, index_col=0)
    dist_matrix = pd.read_csv("./distance_matrices/%s_aa_distance_matrix.csv" % pdb_id, header=0, index_col=0)
    interface = pd.read_csv("./interaction_interface/%s_b_interface.txt" % pdb_id,
                                 names=["chain", "seq", "ins", "name", "interface"])

    interface_values = {"%s_%s_%s" % (row["chain"], row["seq"], row["ins"]): row["interface"] for _, row in interface.iterrows()}
    adj_matrix = graph_weights(dist_matrix, seq_matrix, 4.5)

    print(adj_matrix)

    graph = nx.from_pandas_adjacency(adj_matrix)
    nx.set_node_attributes(G=graph, name="interface", values=interface_values)
    print(nx.get_node_attributes(graph, name="interface"))

    for feature in feature_list:
        feature_values = {"%s_%s_%s" % (row["chain"], row["seq"], row["ins"]):
                              aaindex_features.loc[feature, row["name"]]
                          for _, row in interface.iterrows()}
        nx.set_node_attributes(G=graph, name=feature, values=feature_values)
        print(nx.get_node_attributes(graph, name=feature))

    print(graph.edges.data())
    print(nx.info(graph))

    print(graph.nodes.data())

    nx.draw(graph, with_labels=True, cmap=plt.cm.cool, node_color=[v for v in interface["interface"]])
    plt.show()
    nx.write_gpickle(graph, "%s/%s.gpickle" % (destination_path, pdb_id))
    # graph = nx.read_gpickle("%s.gpickle" % pdb_id)

