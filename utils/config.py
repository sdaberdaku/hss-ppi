import argparse

args = argparse.ArgumentParser()
args.add_argument('--learning_rate', type=float, default=0.001)


args.add_argument('--hidden', type=int, nargs='+', default=[35, 32])

#args.add_argument('--hidden', default=[35, 30,20])
#args.add_argument('--hidden', default=[35, 32])
args.add_argument('--epochs', type=int, default=500)

#args.add_argument('--hidden', default=[38, 38, 38, 38, 18, 18, 18, 18, 9, 9, 6, 6])

args.add_argument('--dropout', type=float, default=0.5)
args.add_argument('--weight_decay', type=float, default=5e-4)
# args.add_argument('--early_stopping', default=10)
# args.add_argument('--max_degree', default=3)

args.add_argument('--protein_class', default='A', choices=["A", "AB", "EI", "ER", "ES", "OG", "OR", "OX"])
args.add_argument('--is_bound', action="store_true", help="If set, consider bound PDB structures, otherwise the unbound ones.")
args.add_argument('--is_receptor', action="store_true", help="If set, consider receptor PDB structures, otherwise the ligand ones.")
args.add_argument('--distance_threshold', type=float, default=12) # 6, 8, 10, 12

args.add_argument('--load_pickle', action="store_true")

args.add_argument('--representation_type', choices=["sequence", "contact_map", "hierarchical"], default="hierarchical", help="Protein representation type.")


args = args.parse_args()
print(args)
