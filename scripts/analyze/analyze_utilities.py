import scripts.parse.parse_utilities
import os
import re

def get_pkl_files(arg):
    """Get file if it is a pkl file."""

    if arg[-1] == "/" or arg[-1] == "\\":
        return [arg + f for f in os.listdir(arg) if scripts.parse.parse_utilities.is_pkl(f)]
    else:
        assert scripts.parse.parse_utilities.is_pkl(arg)
        return [arg]

def get_json_files(arg):
    """Get file if it is a json file."""

    if arg[-1] == "/" or arg[-1] == "\\":
        return [arg + f for f in os.listdir(arg) if scripts.parse.parse_utilities.is_json(f)]
    else:
        assert scripts.parse.parse_utilities.is_json(arg)
        return [arg]

def get_tumor_id(file):
    """Collect tumor ID based on file name."""

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    tumorid = fileName.replace('.pkl', '')
    tumorid = tumorid.replace('.LYSIS.json', '')

    return tumorid

def collect_sumulation_info(simDict, TUMORID, SEED):
    """Populate dictionary with simulation setup information."""

    # Add simulation information to dict
    simDict['TUMOR ID'] = TUMORID
    simDict['SEED'] = int(SEED)

    tsplit = TUMORID.split('_')

    simDict['PLATE'] = tsplit[1]
    simDict['NUTRIENTS'] = 'CONSTANT' if tsplit[0] == 'VITRO' else 'GRAPH VASCULATURE'
    simDict['DOSE'] = int(tsplit[4])
    simDict['TREAT RATIO'] = tsplit[5].replace('-', ':') if tsplit[5] != "NA" else tsplit[5]
    simDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
    simDict['ANTIGENS CANCER'] = int(tsplit[7])
    simDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]

    return simDict

def define_cell_pop_numbers():
    """Define cell population numbers used in simualtions for parsing pkl data."""

    # CELL POP NUMBERS
    POP_CANCER = 0
    POP_HEALTHY = 1
    POP_CD4 = 2
    POP_CD8 = 3

    return POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8

def define_cell_index_numbers():
    """Define cell list index numbers for counting cell population information."""

    # These don't change even if pop number does
    INDEX_CANCER = 0
    INDEX_CANCERLIVE = 1
    INDEX_HEALTHY = 2
    INDEX_HEALTHYLIVE = 3
    INDEX_TCELL = 4
    INDEX_TCELLLIVE = 5
    INDEX_CD4 = 6
    INDEX_CD4LIVE = 7
    INDEX_CD8 = 8
    INDEX_CD8LIVE = 9

    return INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE

def define_pop_names_list():
    """List of all cell population names (where LIVE indicates cells of that type in living cell states)."""

    POP_NAMES = ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE', 'T-CELL', 'T-CELL LIVE', 'CD4', 'CD4 LIVE', 'CD8',
                 'CD8 LIVE']

    return POP_NAMES

def define_pop_indices_list():
    """Define list of cell index numbers."""

    INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE = define_cell_index_numbers()

    POP_INDICES = [INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE,
                   INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE]

    return POP_INDICES

def define_state_names_lists():
    """Define list of possible cell states for tissue cells, living tissue cells, and T-cells."""

    STATES_TISSUE = ['NEUTR', 'APOPT', 'QUIES', 'MIGRA', 'PROLI', 'SENES', 'NECRO']
    STATES_TISSUELIVE = ['NEUTR', 'QUIES', 'MIGRA', 'PROLI', 'SENES']
    STATES_TCELL = ['NEUTR', 'APOPT', 'MIGRA', 'PROLI', 'SENES', 'CYTOT', 'STIMU', 'EXHAU', 'ANERG', 'STARV', 'PAUSE']

    return STATES_TISSUE, STATES_TISSUELIVE, STATES_TCELL

def define_state_indices_lists():
    """Initialize state indicies lists for tissue cells, living tissue cells, and T-cells."""

    STATES_TISSUE_INDICES = [0, 1, 2, 3, 4, 5, 6]

    STATES_TISSUELIVE_INDICES = [0, 2, 3, 4, 5]

    STATES_TCELL_INDICES = [0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 12]

    return STATES_TISSUE_INDICES, STATES_TISSUELIVE_INDICES, STATES_TCELL_INDICES