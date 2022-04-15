from scripts.parse.parse import get_hex_rings
from scripts.parse.parse import load as ABM_load
from scripts.analyze.analyze_utilities import collect_sumulation_info
from scripts.analyze.analyze_utilities import get_pkl_files
from scripts.analyze.analyze_utilities import get_tumor_id
import pickle
import pandas as pd

'''
analyze_env takes a directory of (or a single) .pkl simulation files and
extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:

    TIME
    RADIUS
    SEED
    PLATE
    NUTRIENTS
    DOSE
    TREAT RATIO
    CAR AFFINITY
    ANTIGENS CANCER
    ANTIGENS HEALTHY
    GLUCOSE
    OXYGEN
    TGFA
    IL-2
    GLUCOSE TOTAL
    OXYGEN TOTAL
    TGFA TOTAL
    IL-2 TOTAL
    GLUCOSE TOTAL CONC
    TGFA TOTAL CONC
    IL-2 TOTAL CONC
    
where GLUCOSE, OXYGEN, TGFA, and IL-2 are in the format of a list (per time point) of a list of concentrations at each radius.
TOTAL columns are the total (absolute) amount of molecule in total area at each time point.
TOTAL CONC columns are the total concentration amount of molecule across the entire simulation at each time point.
'''

def make_env_df():
    """Initialize empty dataframe to contain environment information."""

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
               'GLUCOSE',
               'OXYGEN',
               'TGFA',
               'IL-2',
               'GLUCOSE TOTAL',
               'OXYGEN TOTAL',
               'TGFA TOTAL',
               'IL-2 TOTAL',
               'GLUCOSE TOTAL CONC',
               'TGFA TOTAL CONC',
               'IL-2 TOTAL CONC'
    ]
    envDF = pd.DataFrame(columns=columns)

    return envDF

def make_env_dict():
    """Initialize empty dictionary to contain environment information."""

    envDict = { 'TUMOR ID': None,
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
                'GLUCOSE': [],
                'OXYGEN': [],
                'TGFA': [],
                'IL-2': [],
                'GLUCOSE TOTAL': [],
                'OXYGEN TOTAL': [],
                'TGFA TOTAL': [],
                'IL-2 TOTAL': [],
                'GLUCOSE TOTAL CONC': [],
                'TGFA TOTAL CONC': [],
                'IL-2 TOTAL CONC': []
    }

    return envDict

def get_environment_containters(environments):
    """Get enviornment containers from abm_parse environments container output."""

    glucose = environments['glucose']  # fmol/um^3
    oxygen = environments['oxygen']  # mmHg
    tgfa = environments['tgfa']  # pg/cm^3
    IL2 = environments['IL-2']  # molecules IL-2/cm^3

    return glucose, oxygen, tgfa, IL2

def update_analyze_env_dict_with_values_per_radius(envDict, glucose, oxygen, tgfa, IL2, SEED, T, time, HEIGHT):
    """Update enviornment dictionary with species quantities at a given time point."""

    # Add information to dictionary
    envDict['TIME'].append(T[time])
    envDict['GLUCOSE'].append(glucose[SEED][time][HEIGHT])
    envDict['OXYGEN'].append(oxygen[SEED][time][HEIGHT])
    envDict['TGFA'].append(tgfa[SEED][time][HEIGHT])
    envDict['IL-2'].append(IL2[SEED][time][HEIGHT])

    return envDict

def update_analyze_env_dict_with_totals(envDict, glucTotal, oxyTotal, tgfaTotal, IL2Total):
    """Populate enviornment dictionary with total species quantities at a given radius."""

    envDict['GLUCOSE TOTAL'].append(glucTotal)  # fmol
    envDict['OXYGEN TOTAL'].append(oxyTotal)  # average mmHg
    envDict['TGFA TOTAL'].append(tgfaTotal)  # pg
    envDict['IL-2 TOTAL'].append(IL2Total)  # molecules

    return envDict

def update_analyze_env_dict_with_total_concentrations(envDict, glucTotal, tgfaTotal, IL2Total, volTotal):
    """Populate enviornment dictionary with total species concentrations at a given radius."""

    IL2pg = (1E12 * 15500 * IL2Total) / (6.022E23)

    envDict['GLUCOSE TOTAL CONC'].append(glucTotal / volTotal / 1E-12)  # fmol/ml
    envDict['TGFA TOTAL CONC'].append(tgfaTotal / volTotal / 1E-12)  # pg/ml
    envDict['IL-2 TOTAL CONC'].append(IL2pg / volTotal / 1E-12)  # pg/ml

    return envDict

def initiatlize_total_concentration_and_volume_counters():
    """Initialize total concentration counters for each species and volume counter."""

    glucTotal = 0  # fmol
    oxyTotal = 0  # mmHg
    tgfaTotal = 0  # pg
    IL2Total = 0  # molecules

    volTotal = 0  # um^3

    return glucTotal, oxyTotal, tgfaTotal, IL2Total, volTotal

def sum_concentrations_across_simulation(envDict, glucose, oxygen, tgfa, IL2, SEED, time, R, HEIGHT):
    """Calculate concentration of each species across entire simulation radius by calculating total quantity per radius and summing across radii."""

    # Initialize counters
    glucTotal, oxyTotal, tgfaTotal, IL2Total, volTotal = initiatlize_total_concentration_and_volume_counters()

    HEX_VOL_UM = 6780.97  # um^3
    hexRings = get_hex_rings(R)
    radii = [i for i in range(0, len(hexRings))]
    envDict['RADIUS'] = [r + 1 for r in radii]

    for n in range(0, len(hexRings)):
        numHexes = hexRings[n]
        r = radii[n]
        volTotal += numHexes * HEX_VOL_UM
        glucTotal += numHexes * HEX_VOL_UM * glucose[SEED][time][HEIGHT][r]
        oxyTotal += numHexes * HEX_VOL_UM * oxygen[SEED][time][HEIGHT][r]
        tgfaTotal += numHexes * HEX_VOL_UM * 1E-12 * tgfa[SEED][time][HEIGHT][r]
        IL2Total += numHexes * HEX_VOL_UM * 1E-12 * IL2[SEED][time][HEIGHT][r]

    oxyTotal = oxyTotal / volTotal

    return envDict, glucTotal, oxyTotal, tgfaTotal, IL2Total, volTotal

def analyze_env_simulation(envDF, environments, T, R, TUMORID, SEEDS):
    """Collect environment dynamics information for given set of simulations (all seeds)for all time points."""

    HEIGHT = 0

    for SEED in range(0, SEEDS):
        print("\t\t  >" + TUMORID + "_0" + str(SEED) + "_ENVIRONMENT")

        # Make simulation dict
        envDict = make_env_dict()

        # Add simulation information to dict
        envDict = collect_sumulation_info(envDict, TUMORID, SEED)

        glucose, oxygen, tgfa, IL2 = get_environment_containters(environments)

        # Collect environment concentrations at each time
        for time in range(0, len(T)):

            envDict = update_analyze_env_dict_with_values_per_radius(envDict, glucose, oxygen, tgfa, IL2, SEED, T, time, HEIGHT)

            # Calculate totals
            envDict, glucTotal, oxyTotal, tgfaTotal, IL2Total, volTotal = sum_concentrations_across_simulation(envDict, glucose, oxygen, tgfa, IL2, SEED, time, R, HEIGHT)

            envDict = update_analyze_env_dict_with_totals(envDict, glucTotal, oxyTotal, tgfaTotal, IL2Total)

            envDict = update_analyze_env_dict_with_total_concentrations(envDict, glucTotal, tgfaTotal, IL2Total, volTotal)

        # Add tumor information to full simulation dataframe
        envDF = envDF.append(envDict, ignore_index=True)

    return envDF

def analyze_env_simulations(file, TUMORID):
    """Iterate through all files to collect environment dyanmics information."""

    # Make simulation environment dataframe
    envDF = make_env_df()

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    # Load full pickle file
    parsedFile = pickle.load(open(file, "rb"))
    environments = parsedFile['environments']

    envDF = analyze_env_simulation(envDF, environments, T, R, TUMORID, N)

    return envDF

def analyze_env(files, saveLoc):
    """Iterate through all files to pass into functions that collect environment dynamics information."""

    PKLFILES = get_pkl_files(files)

    for file in PKLFILES:

        TUMORID = get_tumor_id(file)

        print(TUMORID)

        envDF = analyze_env_simulations(file, TUMORID)

        if saveLoc != '':
            with open(saveLoc + TUMORID + '_ENVIRONMENT.pkl', 'wb') as f:
                pickle.dump(envDF, f)

    return