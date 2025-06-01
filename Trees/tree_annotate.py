import pandas as pd

# === CONFIGURABLE INPUT ===
mapping_file = "mapping_file.tsv"  # Replace with your actual mapping file path

# === COLOR DEFINITIONS ===
phylum_colors = {
    "Ascomycota": "#E78AC3",       # pink
    "Basidiomycota": "#FDB462",    # orange
    "Mucoromycota": "#BC80BD",     # purple
    "Chytridiomycota": "#80B1D3",  # blue
    "Blastocladiomycota": "#B3DE69",# light green
    "Zoopagomycota": "#A6CEE3",    # pale blue
    "Cryptomycota": "#CCEBC5",     # light turquoise
    "Microsporidia": "#FB8072"     # coral
}

# === LOAD MAPPING FILE ===
df = pd.read_csv(mapping_file, sep="\t")

# === BUILD DATASET_SYMBOL FILE ===
itol_lines = [
    "DATASET_SYMBOL",
    "SEPARATOR TAB",
    "DATASET_LABEL\tPhylum Symbols",
    "COLOR\t#000000",
    "LEGEND_TITLE\tPhylum",
    "LEGEND_SHAPES\t" + "\t".join(["2"] * len(phylum_colors)),  # Using shape '2' (circle) for all
    "LEGEND_COLORS\t" + "\t".join(phylum_colors.values()),
    "LEGEND_LABELS\t" + "\t".join(phylum_colors.keys()),
    "MAXIMUM_SIZE\t20",
    "DATA"
]

# Add each species and its corresponding phylum color
for _, row in df.iterrows():
    phylum = row["Phylum"]
    node_id = row["ID"]
    if phylum in phylum_colors:
        color = phylum_colors[phylum]
        # Format: ID,shape,size,color,fill,position
        itol_lines.append(f"{node_id}\t2\t10\t{color}\t1\t1")

# === WRITE TO FILE ===
with open("phylum_symbols.txt", "w") as f:
    f.write("\n".join(itol_lines))

print("✅ Wrote 'phylum_symbols.txt' — ready for iTOL upload.")

