
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import itertools
from math import isclose
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from scipy.stats import fisher_exact
from itertools import combinations
from matplotlib.colors import ListedColormap
import seaborn as sns


def PTM_matcher(delta_masses, tolerance=0.025):
    """
    Matches observed mass deltas to potential PTM combinations
    """
    # Define PTM masses (same as in R code)
    ptm_list = {
        Methylation = 14.01565,
    Dimethylation = 28.0313,
    Trimethylation = 42.04695,
    Acetylation = 42.0106,
    Propionylation = 56.0262,
    Butyrylation = 70.0419,
    Succinylation = 100.0160,
    Malonylation = 86.0004,
    Formylation = 27.9949,
    Oxidation = 15.9949,
    Dioxidation = 31.9898,
    Nitrosylation = 29.9979,
    Nitration = 44.9851,
    Hydroxylation = 15.9949,
    Phosphorylation = 79.9663,
    Sulfation = 79.9568,
    DisulfideLoss = -2.0157,
    Carbamidomethyl = 57.0215,
    Carboxymethylation = 58.0055,
    Carbamylation = 43.0058,
    Iodoacetamide = 57.0215,
    Hexose = 162.0528,
    HexNAc = 203.0794,
    Deoxyhexose = 146.0579,
    SialicAcid = 291.0954,
    O_Glycosylation = 203.08,
    N_GlycanCore = 1202.42,
    Glycation = 162.0528,
    Palmitoylation = 238.2297,
    Myristoylation = 210.1984,
    Farnesylation = 204.1878,
    Geranylgeranylation = 272.2504,
    Ubiquitination = 114.0429,
    NEDDylation = 114.0429,
    SUMOylation = 383.2281
    }
    
    ptm_combos = []
    for n in range(1, 4):
        ptm_combos.extend(itertools.combinations(ptm_list.keys(), n))
    
    combo_masses = {
        ' + '.join(combo): sum(ptm_list[ptm] for ptm in combo)
        for combo in ptm_combos
    }
    
    results = []
    
    for delta_mass in delta_masses:
        current_tol = tolerance
        trials = 0
        best_match = None
        
        # Adaptive tolerance - increase if no match found
        while trials < 5:
            trials += 1
            # Find matches within current tolerance
            matches = [
                (name, mass) 
                for name, mass in combo_masses.items()
                if isclose(mass, delta_mass, abs_tol=current_tol)
            ]
            
            if matches:
                # Find closest match among candidates
                best_match = min(
                    matches, 
                    key=lambda x: abs(x[1] - delta_mass)
                )
                break
                
            current_tol *= 2  # Double tolerance if no match
        
        if best_match:
            match_name, match_mass = best_match
            results.append({
                'delta_mass': delta_mass,
                'match': match_name,
                'mass': match_mass,
                'delta': match_mass - delta_mass,
                'used_tolerance': current_tol,
                'trials': trials
            })
        else:
            results.append({
                'delta_mass': delta_mass,
                'match': np.nan,
                'mass': np.nan,
                'delta': np.nan,
                'used_tolerance': current_tol,
                'trials': trials
            })
    
    return pd.DataFrame(results)

def parse_pepxml(file_path, sample_size=25):
    """
    Parse pepXML file and extract PSMs with mass deltas
    
    Args:
        file_path: Path to pepXML file
        sample_size: Number of PSMs to sample per chunk
    
    Returns:
        DataFrame with PSM information and mass deltas
    """
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    # Find all spectrum queries
    spectrum_queries = root.findall('.//{http://regis-web.systemsbiology.net/pepXML}spectrum_query')
    total_psms = len(spectrum_queries)
    
    # Create stratified sample indices
    chunk_size = total_psms // 5
    sample_indices = sorted([
        idx 
        for i in range(5)
        for idx in random.sample(
            range(i * chunk_size, min((i + 1) * chunk_size, total_psms)),
            sample_size
        )
    ])
    
    psms = []
    delta_masses = []
    
    for i, spectrum in enumerate(spectrum_queries):
        if i not in sample_indices:
            continue
            
        precursor_mass = float(spectrum.get('precursor_neutral_mass'))
        hits = spectrum.findall('.//{http://regis-web.systemsbiology.net/pepXML}search_hit')
        
        for hit in hits:
            pep_mass = float(hit.get('calc_neutral_pep_mass'))
            delta = precursor_mass - pep_mass
            delta_masses.append(delta)
            
            psms.append({
                'spectrum': spectrum.get('spectrum'),
                'peptide': hit.get('peptide'),
                'protein': hit.get('protein'),
                'precursor_mass': precursor_mass,
                'pep_mass': pep_mass,
                'delta_mass': delta
            })
    
    return pd.DataFrame(psms), delta_masses, total_psms

def plot_ptm_analysis(final_results):
    """Create all PTM analysis plots from R code"""
    plt.figure(figsize=(15, 12))
    
    # Plot 1: PTM Frequency
    plt.subplot(2, 2, 1)
    all_ptms = ' + '.join(final_results['match'].dropna()).split(' + ')
    all_ptms = [ptm.strip().title() for ptm in all_ptms]
    ptm_freq = pd.Series(all_ptms).value_counts().reset_index()
    ptm_freq.columns = ['PTM', 'Count']
    
    plt.bar(ptm_freq['PTM'], ptm_freq['Count'], color='lightblue')
    plt.title('PTM Frequencies')
    plt.ylabel('Count')
    plt.xticks(rotation=90)
    
    # Plot 2: Histogram of PTM combo masses
    plt.subplot(2, 2, 2)
    ptm_masses = final_results['PTM_combo_mass'].dropna()
    bins = np.arange(ptm_masses.min(), ptm_masses.max() + 15, 15)
    plt.hist(ptm_masses, bins=bins, color='lightgreen')
    plt.title('PTM Combo Mass Frequency')
    plt.xlabel('Mass Range')
    plt.ylabel('Count')
    
    # Plot 3: Tolerance usage frequency
    plt.subplot(2, 2, 3)
    tol_freq = final_results['used_tolerance'].value_counts().sort_index()
    tol_freq.plot(kind='bar', color='salmon')
    plt.title('Tolerance Usage Frequency')
    plt.xlabel('Tolerance')
    plt.ylabel('Count')
    
    # Plot 4: Disulfide loss histogram
    plt.subplot(2, 2, 4)
    disulfide = final_results[final_results['match'] == 'DisulfideLoss']['delta_mass']
    plt.hist(disulfide, bins=30, color='orange')
    plt.axvline(x=-2.0157, color='blue', linestyle='--', linewidth=2)
    plt.title('Delta mass deviation ~ Disulfide loss artefact')
    plt.xlabel('Delta Mass (Da)')
    plt.legend(['True Disulfide loss -2.0157 Da'])
    
    plt.tight_layout()
    plt.show()

def co_occurrence_analysis(final_results):
    """Perform PTM co-occurrence analysis and visualization"""
    # Prepare PTM matrix
    ptm_split = final_results['match'].dropna().str.split(' \\+ ')
    unique_ptms = sorted(set().union(*ptm_split))
    
    # Create binary PTM matrix
    ptm_matrix = pd.DataFrame({
        ptm: ptm_split.apply(lambda x: ptm in x)
        for ptm in unique_ptms
    })
    
    # Calculate Jaccard distance
    ptm_dist = pairwise_distances(ptm_matrix.T, metric='jaccard')
    
    # Perform MDS
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_coords = mds.fit_transform(ptm_dist)
    
    # K-means clustering
    k = 6
    km = KMeans(n_clusters=k, random_state=666)
    kmeans_clusters = km.fit_predict(mds_coords)
    
    # Create colormaps
    colors = ListedColormap(sns.color_palette("Set1", n_colors=len(unique_ptms)))
    cluster_colors = ListedColormap(sns.color_palette("Dark2", n_colors=k))
    
    # Create plots
    plt.figure(figsize=(18, 6))
    
    # Plot 1: Text labels
    plt.subplot(1, 3, 1)
    for i, ptm in enumerate(unique_ptms):
        plt.scatter(mds_coords[i, 0], mds_coords[i, 1], color=colors(i), s=100)
        plt.text(mds_coords[i, 0], mds_coords[i, 1], ptm, fontsize=10)
    plt.title('PTM Co-occurrence')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    
    # Plot 2: Colored dots with legend
    plt.subplot(1, 3, 2)
    for i, ptm in enumerate(unique_ptms):
        plt.scatter(mds_coords[i, 0], mds_coords[i, 1], color=colors(i), s=150, label=ptm)
    plt.title('PTM Co-occurrence')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 3: K-means clusters
    plt.subplot(1, 3, 3)
    for i in range(k):
        cluster_pts = mds_coords[kmeans_clusters == i]
        plt.scatter(cluster_pts[:, 0], cluster_pts[:, 1], 
                   color=cluster_colors(i), s=150, label=f'Cluster {i+1}')
    plt.title('MDS K-means Clusters')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.show()
    
    # Fisher's exact test for co-occurrence
    fisher_results = []
    for i, j in combinations(range(len(unique_ptms)), 2):
        ptm1 = unique_ptms[i]
        ptm2 = unique_ptms[j]
        
        a = (ptm_matrix[ptm1] & ptm_matrix[ptm2]).sum()
        b = (ptm_matrix[ptm1] & ~ptm_matrix[ptm2]).sum()
        c = (~ptm_matrix[ptm1] & ptm_matrix[ptm2]).sum()
        d = (~ptm_matrix[ptm1] & ~ptm_matrix[ptm2]).sum()
        
        contingency = [[a, b], [c, d]]
        oddsratio, pvalue = fisher_exact(contingency)
        
        fisher_results.append({
            'pair': f"{ptm1} & {ptm2}",
            'p_value': pvalue,
            'odds_ratio': oddsratio
        })
    
    fisher_df = pd.DataFrame(fisher_results).sort_values('p_value')
    print("\nTop significant PTM co-occurrences:")
    print(fisher_df.head(50))
    
    return fisher_df


if __name__ == "__main__":
    file_path = "testfile.pepXML"
    
    # Parse XML and get PSMs
    psms_df, delta_masses, total_psms = parse_pepxml(file_path)
    
    # Match PTMs
    ptm_matches = PTM_matcher(delta_masses)
    
    # Combine results
    final_results = pd.concat([
        psms_df,
        ptm_matches[['match', 'mass', 'delta', 'used_tolerance', 'trials']]
    ], axis=1)
    
    # Sort by tolerance
    final_results.sort_values('used_tolerance', inplace=True)
    
    final_results.rename(columns={
        'mass': 'PTM_combo_mass',
        'delta': 'delta_fit',
        'trials': '#trials'
    }, inplace=True)
    
    print(final_results[[
        'peptide', 'protein', 'precursor_mass', 'pep_mass',
        'delta_mass', 'PTM_combo_mass', 'delta_fit',
        'match', 'used_tolerance', '#trials'
    ]].to_string(index=False))
    
    print(f"\nTotal sampled PSMs: {len(final_results)} out of {total_psms}")

    plot_ptm_analysis(final_results)
    fisher_results = co_occurrence_analysis(final_results)
    print("Fisher results: ",fisher_results)
