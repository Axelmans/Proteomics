from pyteomics import mzml, pepxml
import os

from src.domain.analyzer import Analyzer

class Reader:

    def __init__(self):
        # use the initializer to define the directory paths
        self.mzml_directory: str = "mzml/"
        self.pepxml_directory: str = "pepxml/"

    def analyse_files(self):
        # 1. iterate over the output pepxml files
        for file in os.listdir(self.pepxml_directory):
            self.__analyze_pepxml(self.pepxml_directory + file)
        # 2.
        pass

    @staticmethod
    def __analyze_mzml(file: str):
        return

    @staticmethod
    def __analyze_pepxml(file: str):
        with pepxml.read(file) as spectra:
            for spectrum in spectra:
                # Each spectrum will be analyzed in analyzer.py
                Analyzer.analyze_spectrum(spectrum)
        return

    @staticmethod
    def __analyse_pin(file: str):
        return
