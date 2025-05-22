
import xml.etree.ElementTree as ET
from Bio import SeqIO
import random
import re

def read_fasta(filepath):
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(filepath, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequences[record.id] = str(record.seq)
    return sequences

def main():
    pepXML_path = "test.pepXML"
    fasta_path = "test.fasta"
    
    tree = ET.parse(pepXML_path)
    root = tree.getroot()
    fasta_dict = read_fasta(fasta_path)
    
    namespace = {'d1': 'http://regis-web.systemsbiology.net/pepXML'}
    spectrum_queries = root.findall('.//d1:spectrum_query', namespace)
    total_psms = len(spectrum_queries)
    
    chunk_size = max(1, total_psms // 5)
    sample_size = 150
    sample_indices = []
    
    for i in range(3):  
        start_idx = i * chunk_size + 1
        end_idx = min((i + 1) * chunk_size, total_psms)
        if start_idx > end_idx:
            continue
        sample_indices.extend(random.sample(range(start_idx, end_idx + 1), 
                             min(sample_size, end_idx - start_idx + 1)))
    
    sample_indices = sorted(sample_indices)
    print(sample_indices)  # Optional: to see which indices were taken
    
    sampled_queries = [spectrum_queries[i-1] for i in sample_indices]  # -1 for 0-based index
    
    # Process sampled PSMs
    psm_data = []
    for query in sampled_queries:
        hit = query.find('.//d1:search_hit', namespace)
        if hit is not None:
            psm_data.append({
                'peptide': hit.get('peptide'),
                'protein': hit.get('protein')
            })
    
    # Find variant sequences - CORE LOGIC
    results = []
    for psm in psm_data:
        peptide = psm['peptide']
        protein = psm['protein']
        
        # Find matching protein sequence (first match)
        protein_seq = None
        for prot_id, seq in fasta_dict.items():
            if protein in prot_id:
                protein_seq = seq
                break
        
        if protein_seq is None:
            continue
            
        peptide_len = len(peptide)
        variants = []
        
        # Slide window through sequence
        for j in range(len(protein_seq) - peptide_len + 1):
            window = protein_seq[j:j+peptide_len]
            
            # Count differences
            diff_count = sum(1 for a, b in zip(peptide, window) if a != b)
            
            # Restrict to only 1 difference
            if diff_count == 1:
                changes = [i for i, (a, b) in enumerate(zip(peptide, window)) if a != b]
                subst = f"{peptide[changes[0]]}→{window[changes[0]]}"
                
                variants.append({
                    'original': peptide,
                    'variant': window,
                    'position': j + changes[0] + 1,  
                    'substitutions': subst
                })
        
        if variants:
            for var in variants:
                results.append({
                    'protein': protein,
                    **var
                })
    
    if results:
        print("\n==== SEQUENCE VARIANTS FOUND ====\n")
        for result in results:
            print(f"Protein: {result['protein']}")
            print(f"Original: {result['original']}")
            print(f"Variant: {result['variant']}")
            print(f"Position: {result['position']}")
            print(f"Substitution: {result['substitutions']}")
            print("-" * 40)
    else:
        print("\nNo sequence variants found in sampled PSMs.")

if __name__ == "__main__":
    main()
