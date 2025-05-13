import pandas as pd
from pyteomics import pepxml
from scipy.stats import chi2_contingency
import glob
import matplotlib.pyplot as plt

def pepxml_to_df(file_path):
    """Convert a pepXML file to a DataFrame of PSMs with error handling."""
    psms = []
    try:
        with pepxml.read(file_path) as reader:
            for psm in reader:
                try:
                    # Safely extract fields with fallbacks
                    spectrum = psm.get("spectrum", "unknown")
                    peptide = psm.get("peptide_sequence", "")
                    
                    # Handle cases where peptideprophet_result might be missing
                    pep_prob = 1.0  # Default to worst-case PEP
                    if "search_hit" in psm and "peptideprophet_result" in psm["search_hit"]:
                        pep_prob = 1 - psm["search_hit"]["peptideprophet_result"].get("probability", 0)
                    
                    psms.append({
                        "spectrum": spectrum,
                        "peptide": peptide,
                        "PEP": pep_prob,
                        "search_type": "open" if "open" in file_path.lower() else "closed"
                    })
                except Exception as psm_error:
                    print(f"Error processing PSM in {file_path}: {psm_error}")
                    continue
    except Exception as file_error:
        print(f"Error reading file {file_path}: {file_error}")
    return pd.DataFrame(psms)

# Load all pepXML files in a folder
files = glob.glob("proteodb/*.pepXML")  # UPDATE THIS PATH
if not files:
    print("No .pepXML files found! Check your path.")
    exit()

print(f"Found {len(files)} .pepXML files to process...")

# Process files with progress feedback
all_psms = []
for i, f in enumerate(files, 1):
    print(f"Processing file {i}/{len(files)}: {f.split('/')[-1]}")
    df = pepxml_to_df(f)
    if not df.empty:
        all_psms.append(df)

if not all_psms:
    print("No valid PSMs found in any files!")
    exit()

all_psms = pd.concat(all_psms)

# Analysis
high_conf = all_psms[all_psms["PEP"] < 0.01]
if high_conf.empty:
    print("No high-confidence PSMs found (PEP < 0.01)!")
    exit()

counts = high_conf["search_type"].value_counts().reset_index()
counts.columns = ["search_type", "high_conf_psms"]

# Statistical test
contingency_table = [
    counts[counts["search_type"] == "open"]["high_conf_psms"].sum(),
    counts[counts["search_type"] == "closed"]["high_conf_psms"].sum()
]
chi2, p, _, _ = chi2_contingency([contingency_table, [len(all_psms)-contingency_table[0], len(all_psms)-contingency_table[1]]])

# Visualization
plt.figure(figsize=(8,6))
bars = plt.bar(counts["search_type"], counts["high_conf_psms"], 
               color=["#1f77b4", "#ff7f0e"], alpha=0.7)

plt.title(f"High-Confidence PSMs by Search Type\n(p-value: {p:.3e})", pad=20)
plt.ylabel("Number of PSMs (PEP < 0.01)")
plt.xlabel("Search Type")

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height):,}',
             ha='center', va='bottom')

plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("open_vs_closed_results.png", dpi=300)
plt.show()

print("\n=== Results ===")
print(counts)
print(f"\nChi-square p-value: {p:.3e}")
if p < 0.05:
    print("The difference is statistically significant (p < 0.05)")
else:
    print("The difference is not statistically significant (p ≥ 0.05)")
