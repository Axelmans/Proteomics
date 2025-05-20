import os
import time
from pyteomics import pepxml

class Analyser:

    def __init__(self):
        self.pepxml_directory: str = "pepxml/"

    def analyse(self):
        start = time.time()
        for number, file in enumerate(os.listdir(self.pepxml_directory)):
            file_start = time.time()
            print(f"Analysing {file} ({number + 1} / {len(os.listdir(self.pepxml_directory))})")
            self.analyze_pepxml(self.pepxml_directory + file)
            print(f"Done Analysing {file} (duration = {(time.time() - file_start): .2f} seconds)")
        print(f"DONE (duration = {(time.time() - start): .2f} seconds)")

    # Read all spectra in a .pepXML file and analyze each one.
    def analyze_pepxml(self, file: str):
        with pepxml.read(file) as spectra:
            for spectrum in spectra:
                # Each spectrum will be analyzed in matcher.py.
                self.analyse_spectrum(spectrum)
        return

    def analyse_spectrum(self, spectrum):
        return

    def results(self):
        return
