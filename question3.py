import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

pep_file = pd.read_csv("pep.tsv", sep="\t")
substitutions = pep_file[
    (pep_file["PEP"] < 0.01) & 
    (pep_file["Peptide"].str.contains("K->R|E->D"))  
]

sub_counts = substitutions["Population"].value_counts().reset_index()
sub_counts.columns = ["Population", "Substitution Count"]

contingency_table = pd.crosstab(substitutions["Population"], columns="Count")
chi2, p, _, _ = chi2_contingency(contingency_table)

sns.barplot(data=sub_counts, x="Population", y="Substitution Count")
plt.title(f"Amino Acid Substitutions by Population\n(p-value: {p:.3f})")
plt.show()