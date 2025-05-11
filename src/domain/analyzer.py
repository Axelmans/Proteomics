from pyteomics import mzml, pepxml
from src.domain.values import ethnic_group, match_ptm_combo

class Analyzer:

    def __init__(self):
        pass

    def analyze_pepxml(self, file: str):
        with pepxml.read(file) as spectra:
            for spectrum in spectra:
                # Each spectrum will be analyzed in analyzer.py
                self.__analyze_spectrum(spectrum)
        return

    # Returns a spectrum as formatted json object
    def __analyze_spectrum(self, spectrum):
        spectrum_data = {
            'group': ethnic_group(spectrum['spectrum']),
            'individual': spectrum['spectrum'][:2],
            'precursor_neutral_mass': spectrum['precursor_neutral_mass'],
            'search_hits': self.__analyze_search_hits(spectrum)
        }
        print(spectrum_data)

    # Analyze search hits separately
    @staticmethod
    def __analyze_search_hits(spectrum):
        analyzed_search_hits = {'search_hits': []}
        precursor_neutral_mass = spectrum['precursor_neutral_mass']
        for hit in spectrum['search_hit']:
            calc_neutral_pep_mass = hit['calc_neutral_pep_mass']
            hit_data = {
                'calc_neutral_pep_mass': calc_neutral_pep_mass,
                'post-translational-modifications': match_ptm_combo(precursor_neutral_mass, calc_neutral_pep_mass)
            }
            analyzed_search_hits['search_hits'].append(hit_data)
        return analyzed_search_hits
