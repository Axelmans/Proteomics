import enum

class EthnicGroups(enum.Enum):

    # These are individual numbers matching to ethnic groups
    AFRICAN = {'02', '10', '32', '48', '50', '87', '103'}
    ASIAN = {'29', '34', '38', '49', '73', '77', '81'}
    HISPANIC = {'13', '14', '27', '41', '50', '57', '60'}
    NATIVE = {'28'}
    OTHER = {'03', '13', '14', '18', '30', '59'}
    WHITE = {'01', '09', '20', '21', '25', '75', '100'}

def ethnic_group(spectrum_name):
    # first 2 characters of the name refer to the number of an individual
    number = spectrum_name[:2]
    # iterate over all the groups and return the one with which the number matches
    groups = []
    for group in EthnicGroups:
        if number in group.value:
            groups.append(group.name)
    # return UNDEFINED in case the number does not match any group (normally this should never happen)
    return groups if len(groups) > 0 else ['UNDEFINED']

class Analyzer:

    @staticmethod
    def analyze_spectrum(spectrum):
        spectrum_data = {
            'group': ethnic_group(spectrum['spectrum']),
            'individual': spectrum['spectrum'][:2]
        }
        # TODO: Add filtering logic (i.e. which spectra are of interest?).
        print(spectrum_data)