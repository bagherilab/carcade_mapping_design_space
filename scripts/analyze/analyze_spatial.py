from scripts.parse.parse import get_hex_rings
from scripts.parse.parse import get_radius
from scripts.parse.parse import load as ABM_load
import scripts.analyze.analyze_utilities
import pickle
import pandas as pd

def make_spatial_df():
    """Initialize empty dataframe to contain cell spatial information."""

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
               'RADIUS',
               'CANCER',
               'CANCER LIVE',
               'CANCER NORMALIZED',
               'CANCER LIVE NORMALIZED',
               'HEALTHY',
               'HEALTHY LIVE',
               'HEALTHY NORMALIZED',
               'HEALTHY LIVE NORMALIZED',
               'T-CELL',
               'T-CELL LIVE',
               'T-CELL NORMALIZED',
               'T-CELL LIVE NORMALIZED',
               'CD4',
               'CD4 LIVE',
               'CD4 NORMALIZED',
               'CD4 LIVE NORMALIZED',
               'CD8',
               'CD8 LIVE',
               'CD8 NORMALIZED',
               'CD8 LIVE NORMALIZED'
    ]
    spatialsDF = pd.DataFrame(columns=columns)

    return spatialsDF

def make_spatial_dict():
    """Initialize empty dictionary to contain cell spatial information."""

    spatialDict = { 'TUMOR ID': None,
                    'SEED': None,
                    'PLATE': None,
                    'NUTRIENTS': None,
                    'DOSE': None,
                    'TREAT RATIO': None,
                    'CAR AFFINITY': None,
                    'ANTIGENS CANCER': None,
                    'ANTIGENS HEALTHY': None,
                    'TIME': [],
                    'RADIUS': [],
                    'CANCER': [],
                    'CANCER LIVE': [],
                    'CANCER NORMALIZED': [],
                    'CANCER LIVE NORMALIZED': [],
                    'HEALTHY': [],
                    'HEALTHY LIVE': [],
                    'HEALTHY NORMALIZED': [],
                    'HEALTHY LIVE NORMALIZED': [],
                    'T-CELL': [],
                    'T-CELL LIVE': [],
                    'T-CELL NORMALIZED': [],
                    'T-CELL LIVE NORMALIZED': [],
                    'CD4': [],
                    'CD4 LIVE': [],
                    'CD4 NORMALIZED': [],
                    'CD4 LIVE NORMALIZED': [],
                    'CD8': [],
                    'CD8 LIVE': [],
                    'CD8 NORMALIZED': [],
                    'CD8 LIVE NORMALIZED': []
    }

    return spatialDict

def make_empty_counts_per_radius_lists(R):
    """Intitalize empty count and normalized by total locations count lists."""

    counts = [[0 for j in range(0, R)] for i in range(0, 10)]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]
    countsNorm = [[0 for j in range(0, R)] for i in range(0, 10)]   # [CANCER NORMALIZED, CANCER LIVE NORMALIZED, HEALTHY NORMALIZED, HEALTHY LIVE NORMALIZED, T-CELL NORMALIZED,
                                                                    # T-CELL LIVE NORMALIZED, CD4 NORMALIZED, CD4 LIVE NORMALIZED, CD8 NORMALIZED, CD8 LIVE NORMALIZED]
    return counts, countsNorm

def count_cells_at_location(agents, time, H, loc, counts, radius):
    """Add cell count based on radius for all cells within a location."""

    POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8 = scripts.analyze.analyze_utilities.define_cell_pop_numbers()

    INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE = scripts.analyze.analyze_utilities.define_cell_index_numbers()

    for pos in range(0, 54):
        if agents[time][H][loc]['pop'][pos] != -1:

            # Find CANCER cells
            if agents[time][H][loc]['pop'][pos] == POP_CANCER:
                counts[INDEX_CANCER][radius] += 1

                # Find CANCER LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CANCERLIVE][radius] += 1

            # Find HEALTHY cells
            if agents[time][H][loc]['pop'][pos] == POP_HEALTHY:
                counts[INDEX_HEALTHY][radius] += 1

                # Find HEALTHY LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_HEALTHYLIVE][radius] += 1

            # Find CD4 cells
            if agents[time][H][loc]['pop'][pos] == POP_CD4:
                counts[INDEX_CD4][radius] += 1
                counts[INDEX_TCELL][radius] += 1

                # Find CD4 LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CD4LIVE][radius] += 1
                    counts[INDEX_TCELLLIVE][radius] += 1

            # Find CD8 cells
            if agents[time][H][loc]['pop'][pos] == POP_CD8:
                counts[INDEX_CD8][radius] += 1
                counts[INDEX_TCELL][radius] += 1

                # Find CD8 LIVE cells
                if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                    counts[INDEX_CD8LIVE][radius] += 1
                    counts[INDEX_TCELLLIVE][radius] += 1

    return counts

def normalize_counts_at_radius_to_radius_locations(counts, countsNorm, radius, radLocs):
    """Normalize count at each radius to number of locations at given radius."""

    POP_INDICES = scripts.analyze.analyze_utilities.define_pop_indices_list()

    for pop in POP_INDICES:
        countsNorm[pop][radius] = counts[pop][radius] / radLocs

    return countsNorm

def update_spatial_dict_with_counts(spatialDict, counts, countsNorm):
    """Populate dictionary containing all cell spatial information with tallied cell information."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    POP_INDICES = scripts.analyze.analyze_utilities.define_pop_indices_list()

    for p in range(0, len(POP_NAMES)):
        spatialDict[POP_NAMES[p]].append(counts[POP_INDICES[p]])
        spatialDict[POP_NAMES[p] + ' NORMALIZED'].append(countsNorm[POP_INDICES[p]])

    return spatialDict

def analyze_spatial_simulation(spatialsDF, agents, T, R, C, TUMORID, SEED):
    """Collect cell spatial dynamics information for given simulation for all time points and locations."""

    # Make simulation dict
    spatialDict = make_spatial_dict()

    # Add simulation information to dict
    spatialDict = scripts.analyze.analyze_utilities.collect_sumulation_info(spatialDict, TUMORID, SEED)
    spatialDict['RADIUS'] = [r + 1 for r in range(0, R)]

    # Set height
    H = 0

    # Get hex rings
    hexRings = get_hex_rings(R)

    for time in range(0, len(T)):

        # Reset list of counts at each radii for each species
        counts, countsNorm = make_empty_counts_per_radius_lists(R)

        for loc in range(0, len(C)):

            location = C[loc]
            radius = get_radius(location)
            radLocs = hexRings[radius]

            # Update cell counts at each location and add to radius count
            counts = count_cells_at_location(agents, time, H, loc, counts, radius)

            countsNorm = normalize_counts_at_radius_to_radius_locations(counts, countsNorm, radius, radLocs)

        # Add information to dictionary
        spatialDict['TIME'].append(T[time])

        spatialDict = update_spatial_dict_with_counts(spatialDict, counts, countsNorm)

    # Add tumor information to full simulation dataframe
    spatialsDF = spatialsDF.append(spatialDict, ignore_index=True)

    return spatialsDF

def analyze_spatial_simulations(file, TUMORID):
    """Iterate through all seeds for a given simulation setup to collect cell spatial dyanmics information."""

    # Make simulation environment dataframe
    spatialsDF = make_spatial_df()

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    for s in range(N):
        print("\t >" + TUMORID + "_0" + str(s) + "_SPATIAL")
        agents = D['agents'][s]
        spatialsDF = analyze_spatial_simulation(spatialsDF, agents, T, R, C, TUMORID, s)

    return spatialsDF

def analyze_spatial(files, saveLoc):
    """Iterate through all files to collect cell dynamcis information.

    analyze_spatial takes a directory of (or a single) .pkl simulation files and
    extracts the data into a dataframe in the form:

        TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

    and saves it to a .pkl where the DATA are the following list of information (also shown in make_spatial_dict and make_spatial_df):

        TIME
        SEED
        PLATE
        NUTRIENTS
        DOSE
        TREAT RATIO
        CAR AFFINITY
        ANTIGENS CANCER
        ANTIGENS HEALTHY
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

    where each cell population is in the format of a list of a list of counts at each radius.

    Usage:
        analyze_spatial(files, saveLoc)

        files
            Path to .pkl files or directory.
        saveLoc
            Location of where to save file.
        sharedLocs
            Collect CAR T-cell information for only CAR T-cells that share a location with at least one cancer cell (default: False).
    """

    PKLFILES = scripts.analyze.analyze_utilities.get_pkl_files(files)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)

    for file in PKLFILES:

        TUMORID = scripts.analyze.analyze_utilities.get_tumor_id(file)

        print(TUMORID)

        spatialDF = analyze_spatial_simulations(file, TUMORID)

        if saveLoc != '':
            with open(saveLoc + TUMORID + '_SPATIAL.pkl', 'wb') as f:
                pickle.dump(spatialDF, f)

    return