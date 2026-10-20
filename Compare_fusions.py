import pandas as pd
from pathlib import Path

def load_fusion_file(*_fusions_combined.txt):
    """Load a fusion file and return a set of fusion pairs."""
    df = pd.read_csv(*_fusions_combined.txt)
    df.drop_duplicates(inplace=True)

    if 'gene1_ensg:gene2_ensg' not in df.columns:
        raise ValueError(f"Missing expected column in: {*_fusions_combined.txt}")
    
    return set(df['gene1_ensg:gene2_ensg'].dropna())

def save_fusions(fusion_set, output_path):
    """Save a set of fusion strings to a file."""
    with open(output_path, 'w') as f:
        f.write("gene1_ensg:gene2_ensg\n")
        for fusion in sorted(fusion_set):
            f.write(f"{fusion}\n")

def compare_fusions(tissue_prefix):
    AD = f"{AD}_unique.txt"
    AD_control = f"{AD_control}_unique.txt"

    try:
        AD = load_fusion_file(short_file)
        AD_control = load_fusion_file(long_file)

        common = AD & AD_control
        unique_to_AD = AD
        unique_to_AD_control = AD_control

        save_fusions(common, f"{tissue_prefix}_common_fusions.txt")
        save_fusions(unique_to_AD, f"{tissue_prefix}_unique_to_AD.txt")
        save_fusions(unique_to_AD_control, f"{tissue_prefix}_unique_to_AD_control.txt")

        print(f"🧬 {tissue_prefix.upper()}:")
        print(f"✅ Common fusions: {len(common)}")
        print(f"🔹 Unique to {AD}: {len(unique_to_short)}")
        print(f"🔸 Unique to {AD_control}: {len(unique_to_long)}")
        print(f"📄 Output: {tissue_prefix}_common_fusions.txt, _unique_to_short.txt, _unique_to_long.txt\n")
    
    except Exception as e:
        print(f"❌ Error comparing fusions for {tissue_prefix}: {e}")

# === Run comparisons for each tissue type ===
for tissue in ['AD', 'AD_control']:
    compare_fusions(tissue)
