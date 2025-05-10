import enum
import itertools

class EthnicGroups(enum.Enum):

    # These are individual numbers matching to ethnic groups
    AFRICAN = {'02', '10', '32', '48', '50', '87', '103'}
    ASIAN = {'29', '34', '38', '49', '73', '77', '81'}
    HISPANIC = {'13', '14', '27', '41', '50', '57', '60'}
    NATIVE = {'28'}
    OTHER = {'03', '13', '14', '18', '30', '59'}
    WHITE = {'01', '09', '20', '21', '25', '75', '100'}

class PTM(enum.Enum):

    # These are known proteomics that will be tried to match (either one or a combination of multiple)
    ACETYLATION = 42.0106
    METHYLATION = 14.0157
    DIMETHYLATION = 28.0314
    TRIMETHYLATION = 42.0470
    PHOSPHORYLATION = 79.9663
    OXIDATION = 15.9949
    CARBAMIDOMETHYLATION = 57.0215
    UBIQUITINATION = 114.0429

def ethnic_group(spectrum_name):
    # first 2 characters of the name refer to the number of an individual
    number = spectrum_name[:2]
    # iterate over all the groups and return the one with which the number matches
    groups = []
    for group in EthnicGroups:
        if number in group.value:
            groups.append(group.name)
    # return UNDEFINED in case the number does not match any group (normally this should never happen)
    return "/".join(groups) if len(groups) > 0 else 'UNDEFINED'

def generate_ptm_combos():
    # generates all combinations of proteomics so mass differences can be matched
    combos = {}
    ptm_items = list(PTM)
    for r in range(1, len(ptm_items) + 1):
        for combo in itertools.combinations(ptm_items, r):
            combos[tuple(ptm.name for ptm in combo)] = sum(ptm.value for ptm in combo)
    return combos

ptm_combos = generate_ptm_combos()

def match_ptm_combo(mass_deficit: int):
    # iterate over all ptm combinations and return the closest match
    closest_match, closest_match_mass = None, None
    for combination, mass in ptm_combos.items():
        if not closest_match or abs(mass_deficit - mass) < abs(mass_deficit - closest_match_mass):
            closest_match, closest_match_mass = combination, mass
    return closest_match
