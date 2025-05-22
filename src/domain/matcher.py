from pyteomics import pepxml
from src.domain.analyser import Analyser
from src.domain.values import get_individual, get_ethnic_group, match_ptm_combo, SEQUENCES
from src.db.db import add_ptm_match, clear_tables

class Matcher(Analyser):

    def __init__(self):
        super().__init__()
        # NOTE: 2 was used for the report, but 0.8 might be a better value for this!
        self.match_diff_threshold = 2
        # Clear the tables before running a new analysis (comment out to prevent data loss!).
        clear_tables()
        self.analyse()

    # Analyzes spectra according to the research questions' needs.
    def analyze_spectrum(self, spectrum):
        for hit in spectrum['search_hit']:
            # We only want hits with rank 1 and a low enough expect score (5%), ignore all others.
            if hit['search_score']['expect'] <= self.expect_threshold:
                # Decoy matches should be ignored.
                if self.__decoy_match(hit):
                    continue
                # Match the mass difference to PTMs.
                individual = get_individual(spectrum['spectrum'])
                ethnic_groups = get_ethnic_group(spectrum['spectrum']).split("/")
                ptm_combo, match_diff = match_ptm_combo(hit['massdiff'])
                if match_diff > self.match_diff_threshold:
                    continue
                for ptm in ptm_combo:
                    for ethnic_group in ethnic_groups:
                        add_ptm_match(individual, ethnic_group, ptm)
        return

    # A decoy match can be recognized by the "rev_" prefix in proteins.
    @staticmethod
    def __decoy_match(hit):
        return any(prot['protein'].startswith("rev_") for prot in hit['proteins'])
