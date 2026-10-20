#!/usr/bin/env python3

import pandas as pd
import argparse
import glob
import sys
import os
from mygene import MyGeneInfo

# ---------------- Step 1: Argument parsing ----------------
parser = argparse.ArgumentParser(description="Extract ENSG IDs from JAFFA output.")
parser.add_argument("sample", help="Sample/tissue name (prefix of the JAFFA file)")
parser.add_argument("--manual-map", default="manual_gene_map.csv",
                    help="Path to CSV with manual rescue gene mappings (gene_symbol,ensembl_id)")
args = parser.parse_args()

sample_name = args.sample
manual_map_file = args.manual_map

# ---------------- Step 2: Find matching JAFFA file ----------------
matches = glob.glob(f"{sample_name}_jaffa.txt")
if not matches:
    print(f"Error: No file found for sample '{sample_name}' matching '*_jaffa.txt'.", file=sys.stderr)
    sys.exit(1)

input_file = matches[0]
output_file = input_file.replace("_jaffa.txt", "_jaffa_fusions_annotated.tsv")

# ---------------- Step 3: Load JAFFA fusions ----------------
fusions = pd.read_csv(input_file, sep="\t")
fusions.columns = fusions.columns.str.strip().str.lower()

# ---------------- Step 4: Extract unique gene names ----------------
gene_names = pd.unique(fusions[['gene1', 'gene2']].values.ravel())

# ---------------- Step 5: Query MyGene.info ----------------
mg = MyGeneInfo()
results = mg.querymany(
    gene_names.tolist(),
    scopes=["symbol", "alias", "name"],
    fields="ensembl.gene",
    species="human"
)

# ---------------- Step 6: Build automated mapping ----------------
symbol_to_ensg = {}
for res in results:
    symbol = res['query']
    if 'ensembl' in res:
        ensembl = res['ensembl']
        if isinstance(ensembl, list):
            symbol_to_ensg[symbol] = ensembl[0]['gene']
        else:
            symbol_to_ensg[symbol] = ensembl['gene']
    else:
        symbol_to_ensg[symbol] = None

# ---------------- Step 7: Load existing manual map ----------------
if os.path.exists(manual_map_file):
    manual_map_df = pd.read_csv(manual_map_file)
    manual_map_dict = dict(zip(manual_map_df['gene_symbol'], manual_map_df['ensembl_id']))
    symbol_to_ensg.update({k: v for k, v in manual_map_dict.items() if pd.notna(v)})
else:
    manual_map_df = pd.DataFrame(columns=["gene_symbol", "ensembl_id"])

# ---------------- Step 8: Annotate fusion file ----------------
fusions['gene1_ensg'] = fusions['gene1'].map(symbol_to_ensg)
fusions['gene2_ensg'] = fusions['gene2'].map(symbol_to_ensg)

# ---------------- Step 9: Save annotated fusion file ----------------
fusions.to_csv(output_file, sep="\t", index=False)

# ---------------- Step 10: Save full mapping table ----------------
mapping_df = pd.DataFrame([
    {"gene_symbol": gene, "ensembl_id": ensg}
    for gene, ensg in symbol_to_ensg.items()
])
mapping_df.to_csv(f"{sample_name}_gene_symbol_to_ensg_map.csv", index=False)

# ---------------- Step 11: Handle still-unmapped genes ----------------
unmapped = [gene for gene in gene_names if symbol_to_ensg.get(gene) is None]

if unmapped:
    # Query MyGene again, stricter (only "symbol")
    strict_results = mg.querymany(
        unmapped,
        scopes="symbol",
        fields="ensembl.gene",
        species="human"
    )

    new_mappings = []
    still_unmapped = []
    for res in strict_results:
        symbol = res['query']
        if 'ensembl' in res:
            ensembl = res['ensembl']
            if isinstance(ensembl, list):
                ensg_id = ensembl[0]['gene']
            else:
                ensg_id = ensembl['gene']
            new_mappings.append({"gene_symbol": symbol, "ensembl_id": ensg_id})
        else:
            still_unmapped.append(symbol)

    # Append only correct new mappings to manual map file
    if new_mappings:
        new_df = pd.DataFrame(new_mappings)
        manual_map_df = pd.concat([manual_map_df, new_df]).drop_duplicates(subset=["gene_symbol"])
        manual_map_df.to_csv(manual_map_file, index=False)
        print(f"📝 Updated manual map file: {manual_map_file} with {len(new_mappings)} new entries")

    # Save remaining unmapped
    if still_unmapped:
        with open(f"{sample_name}_unmapped_genes.txt", "w") as f:
            f.write("\n".join(sorted(still_unmapped)))
        print(f"⚠️ Still unmapped genes saved to: {sample_name}_unmapped_genes.txt")

print(f"✅ Annotated fusion file written to: {output_file}")
print(f"📄 Gene symbol → ENSG mapping saved to: {sample_name}_gene_symbol_to_ensg_map.csv")
