from scripts.parse.parse import load as ABM_load
import scripts.analyze.analyze_utilities
import pickle
import pandas as pd

'''
analyze_cells takes a directory of (or a single) .pkl simulation files
and extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:
    
    TIME
    CANCER
    CANCER LIVE
    HEALTHY
    HEALTHY LIVE
    T-CELL
    T-CELL LIVE
    CD4
    CD4 LIVE
    CD8
    CD8 LIVE
    CANCER LIVE %
    HEALTHY LIVE %
    NEUTR CANCER
    NEUTR CANCER LIVE
    NEUTR HEALTHY
    NEUTR HEALTHY LIVE
    NEUTR T-CELL
    NEUTR CD4
    NEUTR CD8
    NEUTR CANCER %
    NEUTR CANCER LIVE %
    NEUTR HEALTHY %
    NEUTR HEALTHY LIVE %
    NEUTR T-CELL %
    NEUTR CD4 %
    NEUTR CD8 %
    APOPT CANCER
    APOPT HEALTHY
    APOPT T-CELL
    APOPT CD4
    APOPT CD8
    APOPT CANCER %
    APOPT HEALTHY %
    APOPT T-CELL %
    APOPT CD4 %
    APOPT CD8 %
    QUIES CANCER
    QUIES CANCER LIVE
    QUIES HEALTHY
    QUIES HEALTHY LIVE
    QUIES CANCER %
    QUIES CANCER LIVE %
    QUIES HEALTHY %
    QUIES HEALTHY LIVE %
    MIGRA CANCER
    MIGRA CANCER LIVE
    MIGRA HEALTHY
    MIGRA HEALTHY LIVE
    MIGRA T-CELL
    MIGRA CD4
    MIGRA CD8
    MIGRA CANCER %
    MIGRA CANCER LIVE %
    MIGRA HEALTHY %
    MIGRA HEALTHY LIVE %
    MIGRA T-CELL %
    MIGRA CD4 %
    MIGRA CD8 %
    PROLI CANCER
    PROLI CANCER LIVE
    PROLI HEALTHY
    PROLI HEALTHY LIVE
    PROLI T-CELL
    PROLI CD4
    PROLI CD8
    PROLI CANCER %
    PROLI CANCER LIVE %
    PROLI HEALTHY LIVE %
    PROLI HEALTHY %
    PROLI T-CELL %
    PROLI CD4 %
    PROLI CD8 %
    SENES CANCER
    SENES CANCER LIVE
    SENES HEALTHY
    SENES HEALTHY LIVE
    SENES T-CELL
    SENES CD4
    SENES CD8
    SENES CANCER %
    SENES CANCER LIVE %
    SENES HEALTHY %
    SENES HEALTHY LIVE %
    SENES T-CELL %
    SENES CD4 %
    SENES CD8 %
    NECRO CANCER
    NECRO HEALTHY
    NECRO CANCER %
    NECRO HEALTHY %
    CYTOT T-CELL
    CYTOT CD8
    CYTOT T-CELL %
    CYTOT CD8 %
    STIMU T-CELL
    STIMU CD4
    STIMU T-CELL %
    STIMU CD4 %
    EXHAU T-CELL
    EXHAU CD4
    EXHAU CD8
    EXHAU T-CELL %
    EXHAU CD4 %
    EXHAU CD8 %
    ANERG T-CELL
    ANERG CD4
    ANERG CD8
    ANERG T-CELL %
    ANERG CD4 %
    ANERG CD8 %
    STARV T-CELL %
    STARV CD4 %
    STARV CD8 %
    PAUSE T-CELL %
    PAUSE CD4 %
    PUASE CD8%
    AVG CELL CYCLES CANCER
    AVG CELL CYCLES HEALTHY
    AVG CELL CYCLES T-CELL
    AVG CELL CYCLES CD4
    AVG CELL CYCLES CD8
    CELL VOLUMES CANCER
    CELL VOLUMES HEALTHY
    CELL VOLUMES T-CELL
    CELL VOLUMES CD4
    CELL VOLUMES CD8

each in the format of a list of the value of the specified information at each point in time.
'''

def make_cells_df():
    """Initialize empty dataframe to contain cell information."""

    columns = ['TUMOR ID',
               'SEED',
               'PLATE',
               'NUTRIENTS',
               'DOSE',
               'TREAT RATIO',
               'CAR AFFINITY',
               'ANTIGENS CANCER',
               'ANTIGENS HEALTHY',
               'TIME',
               'CANCER',
               'CANCER LIVE',
               'HEALTHY',
               'HEALTHY LIVE',
               'T-CELL',
               'T-CELL LIVE',
               'CD4',
               'CD4 LIVE',
               'CD8',
               'CD8 LIVE',
               'CANCER LIVE %',
               'HEALTHY LIVE %',
               'NEUTR CANCER',
               'NEUTR CANCER LIVE',
               'NEUTR HEALTHY',
               'NEUTR HEALTHY LIVE',
               'NEUTR T-CELL',
               'NEUTR CD4',
               'NEUTR CD8',
               'NEUTR CANCER %',
               'NEUTR CANCER LIVE %',
               'NEUTR HEALTHY %',
               'NEUTR HEALTHY LIVE %',
               'NEUTR T-CELL %',
               'NEUTR CD4 %',
               'NEUTR CD8 %',
               'APOPT CANCER',
               'APOPT HEALTHY',
               'APOPT T-CELL',
               'APOPT CD4',
               'APOPT CD8',
               'APOPT CANCER %',
               'APOPT HEALTHY %',
               'APOPT T-CELL %',
               'APOPT CD4 %',
               'APOPT CD8 %',
               'QUIES CANCER',
               'QUIES CANCER LIVE',
               'QUIES HEALTHY',
               'QUIES HEALTHY LIVE',
               'QUIES CANCER %',
               'QUIES CANCER LIVE %',
               'QUIES HEALTHY %',
               'QUIES HEALTHY LIVE %',
               'MIGRA CANCER',
               'MIGRA CANCER LIVE',
               'MIGRA HEALTHY',
               'MIGRA HEALTHY LIVE',
               'MIGRA T-CELL',
               'MIGRA CD4',
               'MIGRA CD8',
               'MIGRA CANCER %',
               'MIGRA CANCER LIVE %',
               'MIGRA HEALTHY %',
               'MIGRA HEALTHY LIVE %',
               'MIGRA T-CELL %',
               'MIGRA CD4 %',
               'MIGRA CD8 %',
               'PROLI CANCER',
               'PROLI CANCER LIVE',
               'PROLI HEALTHY',
               'PROLI HEALTHY LIVE',
               'PROLI T-CELL',
               'PROLI CD4',
               'PROLI CD8',
               'PROLI CANCER %',
               'PROLI CANCER LIVE %',
               'PROLI HEALTHY LIVE %',
               'PROLI HEALTHY %',
               'PROLI T-CELL %',
               'PROLI CD4 %',
               'PROLI CD8 %',
               'SENES CANCER',
               'SENES CANCER LIVE',
               'SENES HEALTHY',
               'SENES HEALTHY LIVE',
               'SENES T-CELL',
               'SENES CD4',
               'SENES CD8',
               'SENES CANCER %',
               'SENES CANCER LIVE %',
               'SENES HEALTHY %',
               'SENES HEALTHY LIVE %',
               'SENES T-CELL %',
               'SENES CD4 %',
               'SENES CD8 %',
               'NECRO CANCER',
               'NECRO HEALTHY',
               'NECRO CANCER %',
               'NECRO HEALTHY %',
               'CYTOT T-CELL',
               'CYTOT CD8',
               'CYTOT T-CELL %',
               'CYTOT CD8 %',
               'STIMU T-CELL',
               'STIMU CD4',
               'STIMU T-CELL %',
               'STIMU CD4 %',
               'EXHAU T-CELL',
               'EXHAU CD4',
               'EXHAU CD8',
               'EXHAU T-CELL %',
               'EXHAU CD4 %',
               'EXHAU CD8 %',
               'ANERG T-CELL',
               'ANERG CD4',
               'ANERG CD8',
               'ANERG T-CELL %',
               'ANERG CD4 %',
               'ANERG CD8 %',
               'STARV T-CELL',
               'STARV CD4',
               'STARV CD8',
               'STARV T-CELL %',
               'STARV CD4 %',
               'STARV CD8 %',
               'PAUSE T-CELL',
               'PAUSE CD4',
               'PAUSE CD8',
               'PAUSE T-CELL %',
               'PAUSE CD4 %',
               'PAUSE CD8 %',
               'AVG CELL CYCLES CANCER',
               'AVG CELL CYCLES HEALTHY',
               'AVG CELL CYCLES T-CELL',
               'AVG CELL CYCLES CD4',
               'AVG CELL CYCLES CD8',
               'CELL VOLUMES CANCER',
               'CELL VOLUMES HEALTHY',
               'CELL VOLUMES T-CELL',
               'CELL VOLUMES CD4',
               'CELL VOLUMES CD8',
    ]
    cellsDF = pd.DataFrame(columns=columns)

    return cellsDF

def make_cells_dict():
    """Initialize empty dictionary to contain cell information."""

    cellsDict = {   'TUMOR ID': None,
                    'SEED': None,
                    'PLATE': None,
                    'NUTRIENTS': None,
                    'DOSE': None,
                    'TREAT RATIO': None,
                    'CAR AFFINITY': None,
                    'ANTIGENS CANCER': None,
                    'ANTIGENS HEALTHY': None,
                    'TIME': [],
                    'CANCER': [],
                    'CANCER LIVE': [],
                    'HEALTHY': [],
                    'HEALTHY LIVE': [],
                    'T-CELL': [],
                    'T-CELL LIVE': [],
                    'CD4': [],
                    'CD4 LIVE': [],
                    'CD8': [],
                    'CD8 LIVE': [],
                    'CANCER LIVE %': [],
                    'HEALTHY LIVE %': [],
                    'NEUTR CANCER': [],
                    'NEUTR CANCER LIVE': [],
                    'NEUTR HEALTHY': [],
                    'NEUTR HEALTHY LIVE': [],
                    'NEUTR T-CELL': [],
                    'NEUTR CD4': [],
                    'NEUTR CD8': [],
                    'NEUTR CANCER %': [],
                    'NEUTR CANCER LIVE %': [],
                    'NEUTR HEALTHY %': [],
                    'NEUTR HEALTHY LIVE %': [],
                    'NEUTR T-CELL %': [],
                    'NEUTR CD4 %': [],
                    'NEUTR CD8 %': [],
                    'APOPT CANCER': [],
                    'APOPT HEALTHY': [],
                    'APOPT T-CELL': [],
                    'APOPT CD4': [],
                    'APOPT CD8': [],
                    'APOPT CANCER %': [],
                    'APOPT HEALTHY %': [],
                    'APOPT T-CELL %': [],
                    'APOPT CD4 %': [],
                    'APOPT CD8 %': [],
                    'QUIES CANCER': [],
                    'QUIES CANCER LIVE': [],
                    'QUIES HEALTHY': [],
                    'QUIES HEALTHY LIVE': [],
                    'QUIES CANCER %': [],
                    'QUIES CANCER LIVE %': [],
                    'QUIES HEALTHY %': [],
                    'QUIES HEALTHY LIVE %': [],
                    'MIGRA CANCER': [],
                    'MIGRA CANCER LIVE': [],
                    'MIGRA HEALTHY': [],
                    'MIGRA HEALTHY LIVE': [],
                    'MIGRA T-CELL': [],
                    'MIGRA CD4': [],
                    'MIGRA CD8': [],
                    'MIGRA CANCER %': [],
                    'MIGRA CANCER LIVE %': [],
                    'MIGRA HEALTHY %': [],
                    'MIGRA HEALTHY LIVE %': [],
                    'MIGRA T-CELL %': [],
                    'MIGRA CD4 %': [],
                    'MIGRA CD8 %': [],
                    'PROLI CANCER': [],
                    'PROLI CANCER LIVE': [],
                    'PROLI HEALTHY': [],
                    'PROLI HEALTHY LIVE': [],
                    'PROLI T-CELL': [],
                    'PROLI CD4': [],
                    'PROLI CD8': [],
                    'PROLI CANCER %': [],
                    'PROLI CANCER LIVE %': [],
                    'PROLI HEALTHY %': [],
                    'PROLI HEALTHY LIVE %': [],
                    'PROLI T-CELL %': [],
                    'PROLI CD4 %': [],
                    'PROLI CD8 %': [],
                    'SENES CANCER': [],
                    'SENES CANCER LIVE': [],
                    'SENES HEALTHY': [],
                    'SENES HEALTHY LIVE': [],
                    'SENES T-CELL': [],
                    'SENES CD4': [],
                    'SENES CD8': [],
                    'SENES CANCER %': [],
                    'SENES CANCER LIVE %': [],
                    'SENES HEALTHY %': [],
                    'SENES HEALTHY LIVE %': [],
                    'SENES T-CELL %': [],
                    'SENES CD4 %': [],
                    'SENES CD8 %': [],
                    'NECRO CANCER': [],
                    'NECRO HEALTHY': [],
                    'NECRO CANCER %': [],
                    'NECRO HEALTHY %': [],
                    'CYTOT T-CELL': [],
                    'CYTOT CD8': [],
                    'CYTOT T-CELL %': [],
                    'CYTOT CD8 %': [],
                    'STIMU T-CELL': [],
                    'STIMU CD4': [],
                    'STIMU T-CELL %': [],
                    'STIMU CD4 %': [],
                    'EXHAU T-CELL': [],
                    'EXHAU CD4': [],
                    'EXHAU CD8': [],
                    'EXHAU T-CELL %': [],
                    'EXHAU CD4 %': [],
                    'EXHAU CD8 %': [],
                    'ANERG T-CELL': [],
                    'ANERG CD4': [],
                    'ANERG CD8': [],
                    'ANERG T-CELL %': [],
                    'ANERG CD4 %': [],
                    'ANERG CD8 %': [],
                    'STARV T-CELL': [],
                    'STARV CD4': [],
                    'STARV CD8': [],
                    'STARV T-CELL %': [],
                    'STARV CD4 %': [],
                    'STARV CD8 %': [],
                    'PAUSE T-CELL': [],
                    'PAUSE CD4': [],
                    'PAUSE CD8': [],
                    'PAUSE T-CELL %': [],
                    'PAUSE CD4 %': [],
                    'PAUSE CD8 %': [],
                    'AVG CELL CYCLES CANCER': [],
                    'AVG CELL CYCLES HEALTHY': [],
                    'AVG CELL CYCLES T-CELL': [],
                    'AVG CELL CYCLES CD4': [],
                    'AVG CELL CYCLES CD8': [],
                    'CELL VOLUMES CANCER': [],
                    'CELL VOLUMES HEALTHY': [],
                    'CELL VOLUMES T-CELL': [],
                    'CELL VOLUMES CD4': [],
                    'CELL VOLUMES CD8': [],
    }

    return cellsDict

def make_empty_cell_information_lists():
    """Initialize empty lists to contain cell information to be tallied within a time point."""

    counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]

    types = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CANCER
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CANCER LIVE
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # HEALTHY
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # HEALTHY LIVE
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # T-CELL
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # T-CELL LIVE
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD4
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD4 LIVE
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD8
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # CD8 LIVE
             ]

    typesFrac = []

    cycles = [[], [], [], [], [], [], [], [], [], []]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]
    volumes = [[], [], [], [], [], [], [], [], [], []]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]

    return counts, types, typesFrac, cycles, volumes

def caclulate_pop_state_fractions(counts, typesFrac, types):
    """Calculate fractions of cells in each state for all cell populations."""

    # Calculate state fractions for each pop
    for c in range(0, len(counts)):
        if counts[c] == 0:
            typesFrac.append([0 for x in types[c]])
        else:
            typesFrac.append([x / counts[c] for x in types[c]])

    return typesFrac

def update_analyze_cells_dict(cellsDict, time, counts, types, typesFrac, cycles, volumes):
    """Populate dictionary containing all cell information with tallied cell information."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    POP_INDICES = scripts.analyze.analyze_utilities.define_pop_indices_list()

    STATES_TISSUE, STATES_TISSUELIVE, STATES_TCELL = scripts.analyze.analyze_utilities.define_state_names_lists()
    STATES_TISSUE_INDICES, STATES_TISSUELIVE_INDICES, STATES_TCELL_INDICES = scripts.analyze.analyze_utilities.define_state_indices_lists()

    # Add information to dictionary
    cellsDict['TIME'].append(time)

    for p in range(0, len(POP_NAMES)):
        cellsDict[POP_NAMES[p]].append(counts[POP_INDICES[p]])

        if 'LIVE' not in POP_NAMES[p]:

            cellsDict['AVG CELL CYCLES ' + POP_NAMES[p]].append(cycles[POP_INDICES[p]])
            cellsDict['CELL VOLUMES ' + POP_NAMES[p]].append(volumes[POP_INDICES[p]])

            if POP_NAMES[p] == 'CANCER' or POP_NAMES[p] == 'HEALTHY':
                for s in range(0, len(STATES_TISSUE)):
                    key = STATES_TISSUE[s] + ' ' + POP_NAMES[p]
                    cellsDict[key].append(types[POP_INDICES[p]][STATES_TISSUE_INDICES[s]])
                    cellsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TISSUE_INDICES[s]])

            if POP_NAMES[p] == 'T-CELL' or POP_NAMES[p] == 'CD4' or POP_NAMES[p] == 'CD8':
                for s in range(0, len(STATES_TCELL)):
                    if (POP_NAMES[p] == 'CD4' and STATES_TCELL[s] == 'CYTOT'):
                        continue
                    elif (POP_NAMES[p] == 'CD8' and STATES_TCELL[s] == 'STIMU'):
                        continue
                    else:
                        key = STATES_TCELL[s] + ' ' + POP_NAMES[p]
                        cellsDict[key].append(types[POP_INDICES[p]][STATES_TCELL_INDICES[s]])
                        cellsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TCELL_INDICES[s]])

        else:
            if POP_NAMES[p] == 'CANCER LIVE' or POP_NAMES[p] == 'HEALTHY LIVE':
                if counts[POP_INDICES[p - 1]] == 0:
                    cellsDict[POP_NAMES[p] + ' %'].append(0)
                else:
                    cellsDict[POP_NAMES[p] + ' %'].append(counts[POP_INDICES[p]] / counts[POP_INDICES[p - 1]])

                for s in range(0, len(STATES_TISSUELIVE)):
                    key = STATES_TISSUELIVE[s] + ' ' + POP_NAMES[p]
                    cellsDict[key].append(types[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])
                    cellsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])

    return cellsDict

def collect_cell_information(agents, time, H, loc, counts, types, volumes, cycles):
    """Collect all cell information for all cells within a location."""

    POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8 = scripts.analyze.analyze_utilities.define_cell_pop_numbers()
    INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE = scripts.analyze.analyze_utilities.define_cell_index_numbers()

    for pos in range(0, 54):
        if agents[time][H][loc]['pop'][pos] != -1:

            # Find CANCER cells
            if agents[time][H][loc]['pop'][pos] == POP_CANCER:
                counts[INDEX_CANCER] += 1
                types[INDEX_CANCER][agents[time][H][loc]['type'][pos]] += 1
                volumes[INDEX_CANCER].append(agents[time][H][loc]['volume'][pos])
                if agents[time][H][loc]['cycle'][pos] != -1:
                    cycles[INDEX_CANCER].append(agents[time][H][loc]['cycle'][pos])

                # Find CANCER LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CANCERLIVE] += 1
                    types[INDEX_CANCERLIVE][agents[time][H][loc]['type'][pos]] += 1
                    volumes[INDEX_CANCERLIVE].append(agents[time][H][loc]['volume'][pos])
                    if agents[time][H][loc]['cycle'][pos] != -1:
                        cycles[INDEX_CANCERLIVE].append(agents[time][H][loc]['cycle'][pos])

            # Find HEALTHY cells
            if agents[time][H][loc]['pop'][pos] == POP_HEALTHY:
                counts[INDEX_HEALTHY] += 1
                types[INDEX_HEALTHY][agents[time][H][loc]['type'][pos]] += 1
                volumes[INDEX_HEALTHY].append(agents[time][H][loc]['volume'][pos])
                if agents[time][H][loc]['cycle'][pos] != -1:
                    cycles[INDEX_HEALTHY].append(agents[time][H][loc]['cycle'][pos])

                # Find HEALTHY LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_HEALTHYLIVE] += 1
                    types[INDEX_HEALTHYLIVE][agents[time][H][loc]['type'][pos]] += 1
                    volumes[INDEX_HEALTHYLIVE].append(agents[time][H][loc]['volume'][pos])
                    if agents[time][H][loc]['cycle'][pos] != -1:
                        cycles[INDEX_HEALTHYLIVE].append(agents[time][H][loc]['cycle'][pos])

            # Find CD4 cells
            if agents[time][H][loc]['pop'][pos] == POP_CD4:
                counts[INDEX_CD4] += 1
                counts[INDEX_TCELL] += 1
                types[INDEX_CD4][agents[time][H][loc]['type'][pos]] += 1
                types[INDEX_TCELL][agents[time][H][loc]['type'][pos]] += 1
                volumes[INDEX_CD4].append(agents[time][H][loc]['volume'][pos])
                volumes[INDEX_TCELL].append(agents[time][H][loc]['volume'][pos])
                if agents[time][H][loc]['cycle'][pos] != -1:
                    cycles[INDEX_CD4].append(agents[time][H][loc]['cycle'][pos])
                    cycles[INDEX_TCELL].append(agents[time][H][loc]['cycle'][pos])

                # Find CD4 LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CD4LIVE] += 1
                    counts[INDEX_TCELLLIVE] += 1
                    types[INDEX_CD4LIVE][agents[time][H][loc]['type'][pos]] += 1
                    types[INDEX_TCELLLIVE][agents[time][H][loc]['type'][pos]] += 1
                    volumes[INDEX_CD4LIVE].append(agents[time][H][loc]['volume'][pos])
                    volumes[INDEX_TCELLLIVE].append(agents[time][H][loc]['volume'][pos])
                    if agents[time][H][loc]['cycle'][pos] != -1:
                        cycles[INDEX_CD4LIVE].append(agents[time][H][loc]['cycle'][pos])
                        cycles[INDEX_TCELLLIVE].append(agents[time][H][loc]['cycle'][pos])

            # Find CD8 cells
            if agents[time][H][loc]['pop'][pos] == POP_CD8:
                counts[INDEX_CD8] += 1
                counts[INDEX_TCELL] += 1
                types[INDEX_CD8][agents[time][H][loc]['type'][pos]] += 1
                types[INDEX_TCELL][agents[time][H][loc]['type'][pos]] += 1
                volumes[INDEX_CD8].append(agents[time][H][loc]['volume'][pos])
                volumes[INDEX_TCELL].append(agents[time][H][loc]['volume'][pos])
                if agents[time][H][loc]['cycle'][pos] != -1:
                    cycles[INDEX_CD8].append(agents[time][H][loc]['cycle'][pos])
                    cycles[INDEX_TCELL].append(agents[time][H][loc]['cycle'][pos])

                # Find CD8 LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CD8LIVE] += 1
                    counts[INDEX_TCELLLIVE] += 1
                    types[INDEX_CD8LIVE][agents[time][H][loc]['type'][pos]] += 1
                    types[INDEX_TCELLLIVE][agents[time][H][loc]['type'][pos]] += 1
                    volumes[INDEX_CD8LIVE].append(agents[time][H][loc]['volume'][pos])
                    volumes[INDEX_TCELLLIVE].append(agents[time][H][loc]['volume'][pos])
                    if agents[time][H][loc]['cycle'][pos] != -1:
                        cycles[INDEX_CD8LIVE].append(agents[time][H][loc]['cycle'][pos])
                        cycles[INDEX_TCELLLIVE].append(agents[time][H][loc]['cycle'][pos])

    return counts, types, volumes, cycles

def analyze_cell_simulation(cellsDF, agents, T, C, TUMORID, SEED, sharedLocs):
    """Collect cell dynamics information for given simulation for all time points and locations."""

    # Make simulation dict
    cellsDict = make_cells_dict()

    # Add simulation information to dict
    cellsDict = scripts.analyze.analyze_utilities.collect_sumulation_info(cellsDict, TUMORID, SEED)

    POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8 = scripts.analyze.analyze_utilities.define_cell_pop_numbers()

    # Set height
    H = 0

    for time in range(0, len(T)):

        # Reset counts and types
        counts, types, typesFrac, cycles, volumes = make_empty_cell_information_lists()

        for loc in range(0, len(C)):

            if sharedLocs:
                count_cancer = 0
                for pos in range(0, 54):
                    if agents[time][H][loc]['pop'][pos] == POP_CANCER:
                        count_cancer = 1
                        break

                if count_cancer > 0:
                    counts, types, volumes, cycles = collect_cell_information(agents, time, H, loc, counts, types,
                                                                              volumes, cycles)

            else:
                counts, types, volumes, cycles = collect_cell_information(agents, time, H, loc, counts, types, volumes, cycles)

        # Calculate state fractions for each pop
        typesFrac = caclulate_pop_state_fractions(counts, typesFrac, types)

        # Add information to dictionary
        cellsDict = update_analyze_cells_dict(cellsDict, T[time], counts, types, typesFrac, cycles, volumes)

    # Add tumor information to full simulation dataframe
    cellsDF = cellsDF.append(cellsDict, ignore_index=True)

    return cellsDF

def analyze_cell_simulations(file, TUMORID, file_extension, sharedLocs):
    """Iterate through all seeds for a given simulation setup to collect cell dyanmics information."""

    # Make simulation dataframe
    cellsDF = make_cells_df()

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    for s in range(N):
        print("\t >" + TUMORID + "_0" + str(s) + file_extension)
        agents = D['agents'][s]
        cellsDF = analyze_cell_simulation(cellsDF, agents, T, C, TUMORID, s, sharedLocs)

    return cellsDF

def analyze_cells(files, saveLoc, sharedLocs=False):
    """Iterate through all files to collect cell dynamcis information."""

    PKLFILES = scripts.analyze.analyze_utilities.get_pkl_files(files)

    for file in PKLFILES:

        TUMORID = scripts.analyze.analyze_utilities.get_tumor_id(file)

        print(TUMORID)

        if sharedLocs:
            file_extension = '_SHAREDLOCS'

        else:
            file_extension = '_ANALYZED'

        cellsDF = analyze_cell_simulations(file, TUMORID, file_extension, sharedLocs)

        if saveLoc != '':
            with open(saveLoc + TUMORID + file_extension + '.pkl', 'wb') as f:
                pickle.dump(cellsDF, f)

    return