import pandas as pd

colnames = ["ALA", "ARG", "ASN", "ASP", "CYS",
            "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL"]
aaindex1_df = pd.DataFrame(columns=colnames)

with open("aaindex1", "r") as infile:
    text = infile.read()
    entries = text.split("\n//\n")
    for entry in entries:
        fields = entry.split("\n")
        for f in fields:
            if f.startswith("H"):
                aa_index_id = f.split()[1]
        v2 = fields[-1].split()
        v1 = fields[-2].split()
        if ("NA" in v1) or ("NA" in v2):
            continue
        values = [float(x) for x in v1+v2]
        # print(aa_index_id)
        # print(values)
        new_entry = pd.Series(dict(zip(colnames, values)), name=aa_index_id)
        aaindex1_df = aaindex1_df.append(new_entry)

print(aaindex1_df)

aaindex1_df.to_csv("AAIndexFeatures.csv")