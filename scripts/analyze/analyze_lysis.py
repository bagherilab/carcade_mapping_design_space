import scripts.analyze.analyze_utilities
import scripts.parse.parse_utilities
from scripts.parse.parse import get_radius
import pickle
import pandas as pd
'''
analyze_lysis takes a directory of (or a single) .LYSIS.json simulation files
and extracts the data into a dataframe. The resulting file will contain a data frame in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:
    
    TIME
    RADIUS
    TISSUE SEEDED
    CANCER SEEDED
    HEALTHY SEEDED
    T-CELL SEEDED
    CD4 SEEDED
    CD8 SEEDED
    TIME EXACT                              list of exact times of cell death
    RADIUS LYSED EXACT                      list of exact radius cell killed over time 
    TISSUE LYSED EXACT                      list of exact cells killed over time
    TISSUE LYSED CUMULATIVE EXACT           list of cumulative tissue cells killed at exact times
    CANCER LYSED CUMULATIVE EXACT           list of cumulative cancer cells killed at exact times
    HEALTHY LYSED CUMULATIVE EXACT          list of cumulative healthy cells killed at exact times
    TISSUE LYSED TOTAL                      list of total tissue cells killed at collection time points
    CANCER LYSED TOTAL                      list of total cancer cells killed at collection time points
    HEALTHY LYSED TOTAL                     list of total healthy cells killed at collection time points
    
each in the format of a list of the value of the specified information at each point in time.

Usage:
    python analyze_lysis.py FILES [--saveLoc SAVELOC]

    FILES
        Path to .LYSIS.json or directory
    [--saveLoc SAVELOC]
        Location of where to save file, default won't save files
'''

def make_lysis_df():
    """Initialize empty dataframe to contain lysis information."""

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
               'TISSUE SEEDED',
               'CANCER SEEDED',
               'HEALTHY SEEDED',
               'T-CELL SEEDED',
               'CD4 SEEDED',
               'CD8 SEEDED',
               'TIME EXACT',
               'RADIUS LYSED EXACT',
               'TISSUE LYSED EXACT',
               'TISSUE LYSED CUMULATIVE EXACT',
               'CANCER LYSED CUMULATIVE EXACT',
               'HEALTHY LYSED CUMULATIVE EXACT',
               'TISSUE LYSED TOTAL',
               'CANCER LYSED TOTAL',
               'HEALTHY LYSED TOTAL'
    ]
    lysisDF = pd.DataFrame(columns=columns)

    return lysisDF

def make_lysis_dict():
    """Initialize empty dictionary to contain lysis information."""

    lysisDict = {   'TUMOR ID': None,
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
                    'TISSUE SEEDED': [],
                    'CANCER SEEDED': [],
                    'HEALTHY SEEDED': [],
                    'T-CELL SEEDED': [],
                    'CD4 SEEDED': [],
                    'CD8 SEEDED': [],
                    'TIME EXACT': [],
                    'RADIUS LYSED EXACT': [],
                    'TISSUE LYSED EXACT': [],
                    'TISSUE LYSED CUMULATIVE EXACT': [],
                    'CANCER LYSED CUMULATIVE EXACT': [],
                    'HEALTHY LYSED CUMULATIVE EXACT': [],
                    'TISSUE LYSED TOTAL': [],
                    'CANCER LYSED TOTAL': [],
                    'HEALTHY LYSED TOTAL': []
    }

    return lysisDict

def get_simulation_cell_seed_information(lysisDict, LYSIS_JSON):
    """Collect information about simulation setup regarding quantity of each population seeded."""

    POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8 = scripts.analyze.analyze_utilities.define_cell_pop_numbers()

    lysisDict['RADIUS'] = [r for r in range(1, int(LYSIS_JSON['config']['size']['radius']) + 1)]
    lysisDict['CANCER SEEDED'] = int(LYSIS_JSON['config']['pops'][POP_CANCER][-1])
    lysisDict['HEALTHY SEEDED'] = int(LYSIS_JSON['config']['pops'][POP_HEALTHY][-1])
    lysisDict['TISSUE SEEDED'] = lysisDict['CANCER SEEDED'] + lysisDict['HEALTHY SEEDED']

    for helper in LYSIS_JSON['helpers']:
        if helper['type'] == 'TREAT':
            for pop in helper['pops']:
                if pop[0] == POP_CD4:
                    lysisDict['CD4 SEEDED'] = int(pop[1])
                if pop[0] == POP_CD8:
                    lysisDict['CD8 SEEDED'] = int(pop[1])
            lysisDict['T-CELL SEEDED'] = lysisDict['CD4 SEEDED'] + lysisDict['CD8 SEEDED']

    return lysisDict

def get_simulation_radius_information(lysisDict, LYSIS_JSON):
    """Populaon RADIUS dictionary value with list of values of all radii in the simulation."""

    lysisDict['RADIUS'] = [r for r in range(1, int(LYSIS_JSON['config']['size']['radius']) + 1)]

    return lysisDict

def initiatlize_total_cell_counts_per_timepoint():
    """Initialize total cells per timepoint counters."""

    tissueTotal = 0
    cancerTotal = 0
    healthyTotal = 0

    return tissueTotal, cancerTotal, healthyTotal

def initiatlize_exact_lists():
    """Initialize lists containing exact time of cell death and type of cell killed."""

    timeExact = []
    radiusExact = []
    tissueExact = []
    tissueTotalExact = []
    cancerTotalExact = []
    healthyTotalExact = []

    return timeExact, radiusExact, tissueExact, tissueTotalExact, cancerTotalExact, healthyTotalExact

def update_analyze_lysis_dict_with_exact(lysisDict, timeExact, radiusExact, tissueExact, tissueTotalExact, cancerTotalExact, healthyTotalExact):
    """Populate dictionary containing all lysis information with collected exact cell death information."""

    lysisDict['TIME EXACT'] = timeExact
    lysisDict['RADIUS LYSED EXACT'] = radiusExact
    lysisDict['TISSUE LYSED EXACT'] = tissueExact
    lysisDict['TISSUE LYSED CUMULATIVE EXACT'] = tissueTotalExact
    lysisDict['CANCER LYSED CUMULATIVE EXACT'] = cancerTotalExact
    lysisDict['HEALTHY LYSED CUMULATIVE EXACT'] = healthyTotalExact

    return lysisDict

def update_analyze_lysis_dict(lysisDict, timepoint, tissueTotal, cancerTotal, healthyTotal):
    """Populate dictionary containing total cells killed per simualtion timepoint information."""

    # Add information to dictionary
    lysisDict['TIME'].append(timepoint['time'])
    lysisDict['TISSUE LYSED TOTAL'].append(tissueTotal)
    lysisDict['CANCER LYSED TOTAL'].append(cancerTotal)
    lysisDict['HEALTHY LYSED TOTAL'].append(healthyTotal)

    return lysisDict

def count_lysed_cells_per_timepoint(timepoint, LYSIS_JSON, lysisDict):
    """Count how many tissue, cancer, and healthy cells were lysed at each simulation timepoint in json."""

    POP_CANCER, POP_HEALTHY, POP_CD4, POP_CD8 = scripts.analyze.analyze_utilities.define_cell_pop_numbers()

    tissueTotal, cancerTotal, healthyTotal = initiatlize_total_cell_counts_per_timepoint()

    # If timepoint is final timepoint
    if float(timepoint['time']) == float(LYSIS_JSON['config']['days']):

        timeExact, radiusExact, tissueExact, tissueTotalExact, cancerTotalExact, healthyTotalExact = initiatlize_exact_lists()

        for cell in timepoint['cells']:
            timeExact.append(int(cell[0]))
            coord = cell[1]
            coord = coord[0:-2]
            radius = get_radius(coord)
            radiusExact.append(radius)
            tissueExact.append(cell[2][1])

            if cell[2][1] == POP_CANCER:
                cancerTotal += 1
                tissueTotal += 1
                cancerTotalExact.append(cancerTotal)
                healthyTotalExact.append(healthyTotal)
                tissueTotalExact.append(tissueTotal)

            if cell[2][1] == POP_HEALTHY:
                healthyTotal += 1
                tissueTotal += 1
                cancerTotalExact.append(cancerTotal)
                healthyTotalExact.append(healthyTotal)
                tissueTotalExact.append(tissueTotal)

        lysisDict = update_analyze_lysis_dict_with_exact(lysisDict, timeExact, radiusExact, tissueExact, tissueTotalExact, cancerTotalExact, healthyTotalExact)

    else:
        for cell in timepoint['cells']:
            if cell[2][1] == POP_CANCER:
                cancerTotal += 1
                tissueTotal += 1
            if cell[2][1] == POP_HEALTHY:
                healthyTotal += 1
                tissueTotal += 1

    lysisDict = update_analyze_lysis_dict(lysisDict, timepoint, tissueTotal, cancerTotal, healthyTotal)

    return lysisDict

def lysis_sim(lysisDF, LYSIS_JSON, TUMORID):
    """Intitialize lysis dictionary with specific simulation information and iterate through timepoints to collect lysis information."""

    # Make simulation dict
    lysisDict = make_lysis_dict()

    # Add simulation information to dict
    lysisDict = scripts.analyze.analyze_utilities.collect_sumulation_info(lysisDict, TUMORID, LYSIS_JSON['seed'])
    lysisDict = get_simulation_radius_information(lysisDict, LYSIS_JSON)
    lysisDict = get_simulation_cell_seed_information(lysisDict, LYSIS_JSON)

    for timepoint in LYSIS_JSON['timepoints']:

        lysisDict = count_lysed_cells_per_timepoint(timepoint, LYSIS_JSON, lysisDict)

    # Add tumor information to full simulation dataframe
    lysisDF = lysisDF.append(lysisDict, ignore_index=True)

    return lysisDF

def analyze_lysis_simulation(file, TUMORID):
    """Inititalize lysis dictionary and collect lysis file to being processing."""

    # Make simulation environment dataframe
    lysisDF = make_lysis_df()

    # Load tumor
    LYSIS_JSON = scripts.parse.parse_utilities.load_json(file)

    lysisDF = lysis_sim(lysisDF, LYSIS_JSON, TUMORID)

    return lysisDF

def analyze_lysis(files, saveLoc):
    """Iterate through all files to collect lysis dynamcis information.

    analyze_lysis takes a directory of (or a single) .LYSIS.json simulation files and extracts the data into a dataframe. The resulting file will contain a data frame in the form:

        TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

    and saves it to a .pkl where the DATA are the following list of information (also shown in make_lysis_dict and make_lysis_df):

            TIME
            RADIUS
            TISSUE SEEDED
            CANCER SEEDED
            HEALTHY SEEDED
            T-CELL SEEDED
            CD4 SEEDED
            CD8 SEEDED
            TIME EXACT                              list of exact times of cell death
            RADIUS LYSED EXACT                      list of exact radius cell killed over time
            TISSUE LYSED EXACT                      list of exact cells killed over time
            TISSUE LYSED CUMULATIVE EXACT           list of cumulative tissue cells killed at exact times
            CANCER LYSED CUMULATIVE EXACT           list of cumulative cancer cells killed at exact times
            HEALTHY LYSED CUMULATIVE EXACT          list of cumulative healthy cells killed at exact times
            TISSUE LYSED TOTAL                      list of total tissue cells killed at collection time points
            CANCER LYSED TOTAL                      list of total cancer cells killed at collection time points
            HEALTHY LYSED TOTAL                     list of total healthy cells killed at collection time points

    each in the format of a list of the value of the specified information at each point in time.

    Usage:
        analyze_lysis(files, saveLoc)

        files
            Path to .LYSIS.json or directory.
        saveLoc
            Location of where to save file.
    """

    # Get files
    LYSISFILES = scripts.analyze.analyze_utilities.get_json_files(files)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)

    for file in LYSISFILES:

        TUMORID = scripts.analyze.analyze_utilities.get_tumor_id(file)

        print(TUMORID)

        lysisDF = analyze_lysis_simulation(file, TUMORID)

        if saveLoc != '':
            pickle.dump(lysisDF, open(saveLoc + TUMORID + '_LYSED.pkl', "wb"))

    return