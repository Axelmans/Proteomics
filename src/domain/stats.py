import json
import os
from src.domain.analyser import Analyser

class Stats(Analyser):

    def __init__(self):
        super().__init__()
        # the values of interest
        self.spectra = 0
        self.annotations = 0
        self.annotation_ratio = 0
        self.decoy_annotations = 0
        self.target_annotations = 0
        self.target_annotation_ratio = 0
        self.expect_threshold = 0.05
        self.fdr = 0
        self.fdr_threshold = 0.05
        # run analysis upon initializing (comment out if you only want to run it for the closed search)
        # self.analyse()
        # rerun again for the closed search (comment out if you only want to run it for the open search)
        # self.pepxml_directory = "pepxml_closed/"
        self.analyse()
        self.results()

    def analyse_spectrum(self, spectrum):
        self.spectra += 1
        search_hits = spectrum.get('search_hit', [])
        if len(search_hits) == 0:
            return
        for hit in search_hits:
            if hit['hit_rank'] != 1:
                continue
            if hit['search_score']['expect'] <= self.expect_threshold:
                self.annotations += 1
                if self.__decoy_match(hit):
                    self.decoy_annotations += 1
                else:
                    self.target_annotations += 1

    @staticmethod
    def __decoy_match(hit):
        # Decoys can be recognized by their "rev_" prefix
        return any(prot['protein'].startswith("rev_") for prot in hit['proteins'])

    def __reset(self):
        # Resets all values that are modified during runs
        self.spectra = 0
        self.annotations = 0
        self.annotation_ratio = 0
        self.decoy_annotations = 0
        self.target_annotations = 0
        self.target_annotation_ratio = 0
        self.fdr = 0

    def results(self):
        # Store the results in a json file, print a warning if the FDR is too high!
        self.annotation_ratio = round(self.annotations / self.spectra, 3)
        self.target_annotation_ratio = round(self.target_annotations / self.annotations, 3)
        self.fdr = round(self.decoy_annotations / self.target_annotations, 3)
        results_type = "Closed" if "closed" in self.pepxml_directory else "Open"
        results_json = {
            results_type: {
                "Spectra": self.spectra,
                "Expect Threshold": self.expect_threshold,
                "Annotations": self.annotations,
                "Annotation Ratio": self.annotation_ratio,
                "Decoy Annotations": self.decoy_annotations,
                "Target Annotations": self.target_annotations,
                "False Discovery Rate (FDR)": self.fdr,
                "FDR threshold": self.fdr_threshold,
            }
        }
        if self.fdr > self.fdr_threshold:
            print(f"Warning: High FDR detected (> {self.fdr_threshold})")
            print("Consider lowering the expect score threshold!")
        results_file = "src/results/open_stats.json" if results_type == "Open" else "src/results/closed_stats.json"
        with open(results_file, "w") as f:
            json.dump(results_json, f, indent=4)
        self.__reset()

if __name__ == '__main__':
    stats = Stats()
