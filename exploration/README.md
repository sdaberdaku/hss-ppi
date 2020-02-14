# Exploration (January–February 2020)

Scripts from the first experiments that led to HSS-PPI, written before the code in the rest of this repository. They
are kept as they were, for reference. The final pipeline lives in `preprocessing/`, `utils/` and `train.py`.

## Residue features from AAindex

- `format_aaindex.py` reads the AAindex1 database (`aaindex1`) and writes `AAIndexFeatures.csv`: one row per index,
  one column per amino acid. Indices with missing values are skipped.
- `feature_clusters.py` computes the correlation matrix of the indices and clusters them with affinity propagation on
  the absolute correlations. The exemplar of each cluster is written to `selected_features.txt`, so that strongly
  correlated indices are represented once.
- `selected_features.txt` is the 15-index selection from this first run. The lists used by HSS-PPI
  (`data/selected_features_1.txt` and `data/selected_features_2.txt`) were refined later.

## Structure analysis

These scripts run on the bound structures of the protein–protein Docking Benchmark 5 (`benchmark5/structures/bound`)
and on per-residue tables written by earlier versions of the scripts in `preprocessing/` (`relative_ASA/`,
`residue_contacts/`, `distance_matrices/`, `sequence_matrices/`, `interaction_interface/`).

- `check_asa_contacts_correlation.py` measures, for each receptor and ligand, the correlation between relative solvent
  accessibility and the number of residue contacts, and prints the mean and standard deviation over the benchmark.
- `surface_residues_contacts.py` lists the residues with below-average relative accessibility and at least average contact
  counts.
- `aa_CA_distance_matrix.py` writes the residue distance matrix of each structure using C-alpha atoms (falling back to
  the closest heavy atoms when a residue has no C-alpha).
- `build_networkx.py` builds a residue graph for each structure: residues within 4.5 Å of each other, or adjacent in the
  sequence, are connected, and each node carries the interface label and the selected AAindex features. It is an early version of
  the residue graphs used by HSS-PPI.

The scripts use Python 3 with pandas, NumPy, scikit-learn, Biopython and networkx.
