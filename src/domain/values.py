import enum
import itertools
from curses.ascii import isdigit
from pyteomics import fasta

class EthnicGroups(enum.Enum):

    # These are individual numbers matching to ethnic groups
    AFRICAN = {'02', '10', '32', '48', '50', '87', '103'}
    ASIAN = {'29', '34', '38', '49', '73', '77', '81'}
    HISPANIC = {'13', '14', '27', '41', '50', '57', '60'}
    NATIVE = {'28'}
    OTHER = {'03', '13', '14', '18', '30', '59'}
    WHITE = {'01', '09', '20', '21', '25', '75', '100'}

class PTM(enum.Enum):

    # These are known proteomics that will be tried to match (either one or a combination of multiple).
    # Same ordering of Jan's R-project was used.
    METHYLATION = 14.01565
    DIMETHYLATION = 28.0313
    TRIMETHYLATION = 42.04695
    ACETYLATION = 42.0106
    PROPIONYLATION = 56.0262
    BUTYRYLATION = 70.0419
    SUCCINYLATION = 100.0160
    MALONYLATION = 86.0004
    FORMYLATION = 27.9949
    OXIDATION = 15.9949
    DIOXIDATION = 31.9898
    NITROSYLATION = 29.9979
    NITRATION = 44.9851
    HYDROXYLATION = 15.9949
    PHOSPHORYLATION = 79.9663
    SULFATION = 79.9568
    DISULFIDELOSS = -2.0157
    CARBAMIDOMETHYL = 57.0215
    CARBOXYMETHYLATION = 58.0055
    CARBAMYLATION = 43.0058
    IODOACETAMIDE = 57.0215
    HEXOSE = 162.0528
    HEXNAC = 203.0794
    DEOXYHEXOSE = 146.0579
    SIALICACID = 291.0954
    O_GLYCOSYLATION = 203.08
    N_GLYCANCORE = 1202.42
    GLYCATION = 162.0528
    PALMITOYLATION = 238.2297
    MYRISTOYLATION = 210.1984
    FARNESYLATION = 204.1878
    GERANYLGERANYLATION = 272.2504
    UBIQUITINATION = 114.0429
    NEDDYLATION = 114.0429
    SUMOYLATION = 383.2281

def get_individual(spectrum_name):
    # First 3 or 2 characters of the name refer to the number of the individual.
    if spectrum_name[:3].isdigit():
        return spectrum_name[:3]
    elif spectrum_name[:2].isdigit():
        return spectrum_name[:2]
    # Again, this is a case that shouldn't normally occur at all.
    else:
        raise ValueError(f"Couldn't get individual number out of spectrum name: \"{spectrum_name}\"")

def get_ethnic_group(spectrum_name):
    number = get_individual(spectrum_name)
    # Iterate over all the groups and return the one with which the number matches.
    groups = []
    for group in EthnicGroups:
        if number in group.value:
            groups.append(group.name)
    # Return UNDEFINED in case the number does not match any group (normally this should never happen).
    return "/".join(groups) if len(groups) > 0 else 'UNDEFINED'

def __generate_ptm_combos():
    # Generates all combinations of proteomics so mass differences can be matched.
    combos = {}
    ptm_items = list(PTM)
    # Limit to combinations of at most 3 items.
    for r in range(1, 4):
        for combo in itertools.combinations(ptm_items, r):
            combos[tuple(ptm.name for ptm in combo)] = sum(ptm.value for ptm in combo)
    return combos

ptm_combos = __generate_ptm_combos()

def match_ptm_combo(mass_diff: float):
    # iterate over all ptm combinations and return the closest match.
    best_match = list(min(ptm_combos.items(), key=lambda item: abs(item[1] - mass_diff))[0])
    return best_match, abs(mass_diff - ptm_combos[tuple(best_match)])
