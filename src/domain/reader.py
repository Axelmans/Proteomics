from pyteomics import mzml, pepxml

import matplotlib
import matplotlib.pyplot as plt
import os

matplotlib.use('Agg')


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
        pass

    @staticmethod
    def __analyze_pepxml(file: str):
        with pepxml.read(file) as spectra:
            for spectrum in spectra:
                print(spectrum)
                annotation = spectrum['spectrum']
                search_hit = spectrum['search_hit']

    @staticmethod
    def __analyse_pin(file: str):
        pass
