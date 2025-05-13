import os
from src.domain.analyzer import Analyzer

class Reader:

    def __init__(self):
        # use the initializer to define the directory paths
        self.analyzer = Analyzer()
        self.pepxml_directory: str = "pepxml/"

    def read_files(self):
        # 1. iterate over the output pepxml files
        for file in os.listdir(self.pepxml_directory):
            self.analyzer.analyze_pepxml(self.pepxml_directory + file)
            break
