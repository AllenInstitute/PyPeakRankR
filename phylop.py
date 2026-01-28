import pandas as pd
import pyBigWig

# ---------- PARAMETERS ----------
PHYLOP_BW = "hg38.phyloP100way.bw"
HUMAN_BED = "human_peaks_fixed.bed"
MACAQUE_BED = "peaks_ucsc_all.bed"
OUTPUT = "peaks_NC_with_phyloP.tsv"
MAX_LIFTOVER_LEN = 5000  # anything longer = 0 phyloP

# ---------- LOAD ORIGINAL MACAQUE BED ----------
# Automatically detect all columns
macaque_bed = pd.read_csv(MACAQUE_BED, sep="\t", header=None)
macaque_bed.columns = [f"col{i+1}" for i in range(macaque_bed.shape[1])]

# ---------- LOAD HUMAN LIFTOVER BED ----------
human_bed = pd.read_csv(HUMAN_BED, sep="\t", header=None)
human_bed.columns = ["chr_human", "start_human", "end_human", "group", "status", "name"]

# ---------- MERGE BY PEAK NAME ----------
# Assumes 'name' is always the 6th column in macaque_bed (like AiE1211q)
merged = macaque_bed.merge(
    human_bed[["chr_human", "start_human", "end_human", "name"]],
    left_on="col6", right_on="name", how="left"
)

# ---------- OPEN PHYLOP BIGWIG ----------
bw = pyBigWig.open(PHYLOP_BW)

# ---------- COMPUTE MEAN PHYLOP ----------
phyloP_scores = []
for _, row in merged.iterrows():
    if pd.isna(row["chr_human"]):
        phyloP_scores.append(0)
        continue

    start_h = int(row["start_human"])
    end_h = int(row["end_human"])

    # skip abnormally long liftOver regions
    if end_h - start_h > MAX_LIFTOVER_LEN:
        phyloP_scores.append(0)
        continue

    try:
        vals = bw.values(row["chr_human"], start_h, end_h)
        vals = [v for v in vals if v is not None]
        mean_score = sum(vals) / len(vals) if vals else 0
    except Exception:
        mean_score = 0

    phyloP_scores.append(mean_score)

bw.close()

# ---------- ADD PHYLOP COLUMN ----------
merged["phyloP_mean"] = phyloP_scores

# ---------- CONVERT UCSC CHR TO NC IDs ----------
# Load mapping table: UCSC → NC IDs
mapping = pd.read_csv("rheMac10.chromAlias.txt", sep="\t", comment="#", header=None)
mapping.columns = ["sequenceName", "alias_names", "chr_nc"]
mapping = mapping[["sequenceName", "chr_nc"]]

# Merge to get NC IDs for macaque chr
merged = merged.merge(mapping, left_on="col1", right_on="sequenceName", how="left")
merged["col1"] = merged["chr_nc"]
merged = merged.drop(columns=["sequenceName","chr_nc"])

# ---------- OUTPUT ----------
# Keep *exactly* the original columns + phyloP_mean at the end
output_cols = list(macaque_bed.columns) + ["phyloP_mean"]
output = merged[output_cols]

output.to_csv(OUTPUT, sep="\t", index=False, header=False)
