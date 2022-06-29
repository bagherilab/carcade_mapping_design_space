from scripts.analyze.analyze_cells import make_cells_df
from scripts.analyze.analyze_env import make_env_df
from scripts.analyze.analyze_spatial import make_spatial_df
from scripts.analyze.analyze_lysis import make_lysis_df

def make_datatype_specific_df(TYPE):
    """Make dataframes for treated and untreated data based on file type."""

    # Make dataframe based on type given
    if TYPE == 'ENVIRONMENT':
        simsDF = make_env_df()
        untreatedDF = make_env_df()
    elif TYPE == 'SPATIAL':
        simsDF = make_spatial_df()
        untreatedDF = make_spatial_df()
    elif TYPE == 'LYSED':
        simsDF = make_lysis_df()
        untreatedDF = make_lysis_df()
    elif TYPE == 'SHAREDLOCS':
        simsDF = make_cells_df()
        untreatedDF = make_cells_df()
    else:
        if TYPE != 'ANALYZED':
            print('No valid type specified. Assuming ANALYZE.')
            TYPE = 'ANALYZED'
        simsDF = make_cells_df()
        untreatedDF = make_cells_df()

    return simsDF, untreatedDF, TYPE

def make_list_columns_list():
    """Make list of columns in dataframe with list of list-values."""

    LIST_COLUMNS = ["AVG CELL CYCLES CANCER", "AVG CELL CYCLES HEALTHY", "AVG CELL CYCLES T-CELL",
                    "AVG CELL CYCLES CD4", "AVG CELL CYCLES CD8",
                    "CELL VOLUMES CANCER", "CELL VOLUMES HEALTHY", "CELL VOLUMES T-CELL", "CELL VOLUMES CD4",
                    "CELL VOLUMES CD8"]

    return LIST_COLUMNS

def make_options_dict():
    """Inititlaize empty options dictionary to help name file based on subset requested where X indicates all values of that feature present in subset."""

    # Set up options dictionary for naming save file
    optionsDict = {
        'DOSE': 'X',
        'TREAT RATIO': 'X',
        'CAR AFFINITY': 'X',
        'ANTIGENS CANCER': 'X',
        'ANTIGENS HEALTHY': 'X'
    }

    return optionsDict

def construct_file_save_name(optionsDict):
    """Construct file save name based on options dictionary."""

    # Construct save file name
    name = optionsDict['DOSE'] + '_' + optionsDict['TREAT RATIO'] + '_' + \
           optionsDict['CAR AFFINITY'] + '_' + optionsDict['ANTIGENS CANCER'] + '_' + \
           optionsDict['ANTIGENS HEALTHY'] + '_'

    return name

def print_save_message(subsetsRequested):
    """Print final save message after data is successfully saved."""

    if subsetsRequested == '':
        print('Data saved!')
    else:
        print('Set saved!')

    return