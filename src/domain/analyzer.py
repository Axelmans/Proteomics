from src.domain.values import ethnic_group, match_ptm_combo

class Analyzer:

    @staticmethod
    def analyze_spectrum(spectrum):
        spectrum_data = {
            'group': ethnic_group(spectrum['spectrum']),
            'individual': spectrum['spectrum'][:2],
            'precursor_neutral_mass': spectrum['precursor_neutral_mass'],
            'search_hit': spectrum['search_hit']
        }
        print(spectrum_data)