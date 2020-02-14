import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation

features = pd.read_csv("AAIndexFeatures.csv", header=0, index_col=0)
features = features.transpose()
feature_names = features.columns.values
correlation_matrix = features.corr()
correlation_matrix.to_csv("feature_correlation_matrix.csv")

clustering = AffinityPropagation(affinity="precomputed", verbose=True, preference=-2).fit((np.abs(correlation_matrix.values) - 1))
labels = clustering.labels_
cluster_centers = clustering.cluster_centers_indices_
print(cluster_centers)
with open("selected_features.txt", "w") as features_list:
    for i in range(0, len(cluster_centers)):
        features_list.write("%s\n" % feature_names[clustering.cluster_centers_indices_[i]])