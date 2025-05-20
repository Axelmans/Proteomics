from pyteomics import pepxml
from src.domain.analyser import Analyser
from src.domain.values import get_individual, get_ethnic_group, match_ptm_combo, SEQUENCES
from src.db.db import add_ptm_match, clear_tables

class Matcher(Analyser):

    def __init__(self):
        super().__init__()
        # Clear the tables before running a new analysis (commented out to prevent data loss).
        # clear_tables()
        self.analyse()

    # Read all spectra in a .pepXML file and analyze each one.
    def analyze_pepxml(self, file: str):
        with pepxml.read(file) as spectra:
            for spectrum in spectra:
                # Each spectrum will be analyzed in matcher.py.
                self.analyze_spectrum(spectrum)
        return

    # Analyzes spectra according to the research questions' needs.
    def analyze_spectrum(self, spectrum):
        self.spectra += 1
        if len(spectrum['search_hit']) == 0:
            print("No hits")
        else:
            for hit in spectrum['search_hit']:
                # We only want hits with rank 1 and a low enough expect score (5%), ignore all others.
                if hit['search_score']['expect'] <= self.expect_threshold:
                    # Decoy matches should be ignored.
                    if self.decoy_match(hit):
                        continue
                    # QUESTION 2: Match mass differences to PTMs.
                    individual = get_individual(spectrum['spectrum'])
                    ethnic_groups = get_ethnic_group(spectrum['spectrum']).split("/")
                    ptm_combo, match_diff = match_ptm_combo(hit['massdiff'])
                    if match_diff > 2:
                        continue
                    for ptm in ptm_combo:
                        for ethnic_group in ethnic_groups:
                            add_ptm_match(individual, ethnic_group, ptm)
                    self.annotated_spectra += 1
        return

    # Prints the stats for the first research question
    def get_annotation_stats(self):
        print("Spectra Annotations:")
        print(f"Amount: {self.annotated_spectra}")
        if self.spectra == 0:
            print("Percentage: 0.00% (No spectra analysed)")
        else:
            percentage = self.annotated_spectra / self.spectra * 100
            print(f"Percentage: {percentage:.2f}")

    # A decoy match can be recognized by the "rev_" prefix in proteins.
    def decoy_match(self, hit):
        if any(prot['protein'].startswith("rev_") for prot in hit['proteins']):
            self.decoy_match_count += 1
            return True
        else:
            self.target_match_count += 1
            return False
