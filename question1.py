import pandas as pd
from pyteomics import pepxml
from scipy.stats import chi2_contingency
import glob
import matplotlib.pyplot as plt

def pepxml_to_df(file_path):
    """Convert a pepXML file to a DataFrame of PSMs."""
    psms = []
    with pepxml.read(file_path) as reader:
        for psm in reader:
            psms.append({
                "spectrum": psm["spectrum"],
                "peptide": psm["peptide_sequence"],
                "PEP": 1 - psm["search_hit"]["peptideprophet_result"]["probability"],  # Convert to PEP
                "search_type": "open" if "open" in file_path.lower() else "closed"  # Classify search type
            })
    return pd.DataFrame(psms)

files = glob.glob("LocalAll/*.pepXML")  
all_psms = pd.concat([pepxml_to_df(f) for f in files])
high_conf = all_psms[all_psms["PEP"] < 0.01]
counts = high_conf["search_type"].value_counts().reset_index()
counts.columns = ["search_type", "high_conf_psms"]
print(counts)

contingency_table = [
    [len(high_conf[high_conf["search_type"] == "open"]), 
     len(high_conf[high_conf["search_type"] == "closed"])],
    [len(all_psms) - len(high_conf), len(all_psms)]  # Low-confidence PSMs 
]
chi2, p, _, _ = chi2_contingency(contingency_table)
print(f"Chi-square p-value: {p:.3f}")

plt.bar(counts["search_type"], counts["high_conf_psms"], color=["blue", "orange"])
plt.title(f"High-Confidence PSMs by Search Type\n(p-value: {p:.3f})")
plt.ylabel("Number of PSMs (PEP < 0.01)")
plt.xlabel("Search Type")
plt.show()