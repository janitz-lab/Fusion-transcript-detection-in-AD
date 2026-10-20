###---------------STAR-Fusion----------------------###

import pandas as pd
import re
import argparse
import glob
import sys


# -------------Step 1: Argument parsing----------------------

parser = argparse.ArgumentParser(description="Extract ENSG IDs from STAR-Fusion output.")
parser.add_argument("sample", help="Sample/tissue name (prefix of the STAR-Fusion file)")
args = parser.parse_args()

sample_name = args.sample

# -----------Step 2: Find matching STAR-Fusion file-----------
matches = glob.glob(f"{sampmle_name}_combined_star_fusions.txt")
if not matches:
    print(f"Error: No file found for sample '{sample_name}' matching '*_combined_star_fusions.txt'.", file=sys.stderr)
    sys.exit(1)

input_file = matches[0]
output_file = input_file.replace("_combined_star_fusions.txt", "_combined_star_fusions_annotated.txt")


# -------------Step 3: Load STAR-Fusion output-------------
df = pd.read_csv(input_file, sep="\t", header=None)

# Extract columns containing ENSG IDs
left_col = df[7] # e.g., 'SEPTIN7P14^ENSG00000245958.6'
right_col = df[9] # e.g., 'PSPH^ENSG00000146733.14'

# ---Step 4:Function to extract base ENSG without version-----
def extract_ensg(entry):
    if pd.isna(entry):
        return None
    match = re.search(r"ENSG\d+", entry)
    return match.group(0) if match else None

# --Apply extraction--
df["gene1_ENSG"] = left_col.apply(extract_ensg)
df["gene2_ENSG"] = right_col.apply(extract_ensg)

# -------------Step 5:Save updated dataframe--------------------
df.to_csv(output_file, sep="\t", index=False, header=False)
print(f"Processed: {input_file} → {output_file}")

