#!/usr/bin/env python3
####___________________ARRIBA________________________#####

import pandas as pd
import re
import argparse
import glob
import sys

# ----------------Step 1: Argument parsing---------------

parser = argparse.ArgumentParser(description="Extract ENSG IDs from Arriba output.")
parser.add_argument("sample", help="Sample/tissue name (prefix of the Arriban file)")
args = parser.parse_args()

sample_name = args.sample

# ---------Step 2: Find matching Arriba file---------

matches = glob.glob(f"{sample_name}_arriba_fusions_high.txt")
if not matches:
    print(f"Error: No file found for sample '{sample_name}' matching '*_arriba_fusions_high.txt'.", file=sys.stderr)
    sys.exit(1)

input_file = matches[0]
output_file = input_file.replace("_arriba_fusions_high.txt", "_arriba_fusions_annotated.txt")

#--------------Step 3: Load Arriba output------------------
arriba_df = pd.read_csv(input_file, sep='\t', header=None)

# -----Strip ENSG version (keep only before the dot)-----
arriba_df_cleaned = arriba_df.apply(
    lambda col: col.map(
        lambda x: x.split('.')[0] if isinstance(x, str) and x.startswith("ENSG") else x
    )
)

#--------------Step 4:Save cleaned file--------------------
arriba_df_cleaned.to_csv(output_file, sep='\t', index=False, header=False)

print(f"Cleaned Arriba file saved as '{output_file}'")



