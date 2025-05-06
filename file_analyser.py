from pyteomics import mzml, pepxml

import matplotlib
import matplotlib.pyplot as plt
import os

matplotlib.use('Agg')


class FileAnalyser:

    def __init__(self):
        # use the initializer to define the directory paths
        self.directory: str = os.path.dirname(os.path.abspath(__file__))
        self.graph_directory: str = self.directory + "/graphs/"
        self.input_directory: str = self.directory + "/input/"
        self.results_directory: str = self.directory + "/results/"

    def analyse_files(self):
        # 1. iterate over the input files in the "input" directory
        for file in os.listdir(self.input_directory):
            full_path = os.path.join(self.input_directory, file)
            # .mzml files
            if os.path.isfile(full_path) and full_path.lower().endswith(".mzml"):
                self.__analyze_mzml(full_path)
        # 2. iterate over all the files in the "results" directory
        for file in os.listdir(self.results_directory):
            full_path = os.path.join(self.results_directory, file)
            # .pepxml files
            if os.path.isfile(full_path) and full_path.lower().endswith(".pepxml"):
                self.__analyze_pepxml(full_path)
            # .pin files
            pass

    def __analyze_mzml(self, file: str):
        with mzml.read(file) as spectra:
            for i, spectrum in enumerate(spectra):
                # Fetch info from the file
                # Reference: https://pyteomics.readthedocs.io/en/latest/api/mzml.html
                mz = spectrum['m/z array']
                intensity = spectrum.get('intensity array', [])
                scan_id = spectrum.get('id', f"spectrum_{i}")
                # Create a plot using matplotlib and save it as a .png file inside the "graphs" directory
                plt.figure(figsize=(10, 4))
                plt.plot(mz, intensity, color='blue', linewidth=1)
                plt.title(f"Spectrum: {scan_id}")
                plt.xlabel("m/z")
                plt.ylabel("Intensity")
                plt.tight_layout()
                plt.savefig(file + ".png")
                break

    @staticmethod
    def __analyze_pepxml(file: str):
        pass

    @staticmethod
    def __analyse_pin(file: str):
        pass
