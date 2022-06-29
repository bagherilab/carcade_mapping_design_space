import scripts.analyze.analyze_utilities
import scripts.subset.subset_utilities
import os
import pickle

__author__ = "Alexis N. Prybutok"
__email__ = "aprybutok@u.northwestern.edu"

'''
ABM_SUBSET takes a directory of (or a single) .pkl analyzed files
and extracts the specified subset of the data into a dataframe in the form of the
corresponding data file type (analyze, lysis, environment, spatial).

Usage:
    python subset.py FILES XMLNAME TYPE SAVELOC SUBSET STATES

    FILES
        Path to .pkl or directory
    XMLNAME
        Name of XML file used to generate dataset for purpose of saving
    TYPE
        Type of analysis file given (analyzed, environment, spatial, or lysis)
    SAVELOC
        Location of where to save file, default will save here
    SUBSET
        Feature(s) by which to subselect data and plot only subset of data
        with given specified value (ex: [DOSE,500] will only plot all data
        with a dose of 500 cells) given as semicolon separated list of sets of 
        feature name and value of feature by which to plot separated by a colon.
        More than one feature can be chosen per plot, as indicated by commas 
        between sets, but each new list will make a different set of plots. (ex: 
        [(DOSE:500),(TREAT RATIO:50-50)];[(DOSE:500),(CAR ANTIGEN:1000)] will 
        make two sets of plot, where the first is only data with both DOSE 
        500 and TREAT RATIO 50:50, and the second is only data with DOSE 500
        and ANTIGEN 1000) (default: none, all data plotted on one graph)
    STATES
        Flag indicating only to save state data (not list column data) (default: False).
'''
def parse_requested_subsets(subsetsRequested):
    """Parse subsets requested based on subset input."""

    # Parse feature plotting
    SUBSETS = []
    if subsetsRequested != '':
        ssplit = str(subsetsRequested).split(";")
        for s in ssplit:
            subset = []
            psplit = s.split(",")
            for p in psplit:
                vsplit = p.split(":")
                label = vsplit[0].strip('[]()')
                value = vsplit[1].strip('[]()')
                if label == 'ANTIGENS CANCER' or label == 'DOSE':
                    value = int(value)
                elif label == 'CAR AFFINITY':
                    value = float(value)
                elif label == 'TREAT RATIO':
                    value = value.replace('-', ':')
                subset.append([label, value])
            SUBSETS.append(subset)

    return SUBSETS

def check_if_file_in_subset(file, files, optionsDict):
    """Determine if a file should be included in a subset based on file feature values."""

    file_in_subset = False

    filename = file.replace(files, '')
    fileSplit = filename.split("_")

    # Check file name for set specifications and open file and add to set only if match
    if fileSplit[4] == optionsDict['DOSE'] or optionsDict['DOSE'] == 'X':
        if fileSplit[5] == optionsDict['TREAT RATIO'] or optionsDict['TREAT RATIO'] == 'X':
            try:
                optionsAffinity = float(optionsDict['CAR AFFINITY'])
            except ValueError:
                optionsAffinity = 'X'
            if float(fileSplit[6]) == optionsAffinity or optionsDict['CAR AFFINITY'] == 'X':
                if fileSplit[7] == optionsDict['ANTIGENS CANCER'] or optionsDict['ANTIGENS CANCER'] == 'X':
                    if fileSplit[8] == optionsDict['ANTIGENS HEALTHY'] or optionsDict['ANTIGENS HEALTHY'] == 'X' or fileSplit[8] == 'NA':
                        file_in_subset = True
                        print('\t' + filename)

    return file_in_subset

def drop_list_columns(simsDF, untreatedDF):
    """Drop columns in dataframe where values are lists (ex: volume, cell cycle legth)."""

    LIST_COLUMNS = scripts.subset.subset_utilities.make_list_columns_list()

    simsDF = simsDF.drop(LIST_COLUMNS, axis=1)
    untreatedDF = untreatedDF.drop(LIST_COLUMNS, axis=1)

    return simsDF, untreatedDF

def retrieve_file_data(file, states):
    """Load file information from given pkl based on if only states (and not list columns) selected."""

    with open(file, 'rb') as f:
        simDF = pickle.load(f)
    if states:
        simDF = simDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']

    return simDF

def update_options_dict_antigens_healthy_if_no_healthy_cells(optionsDict, untreatedName, TYPE):
    """Populate options dictionary with NA for healthy cell antigen value if healthy cells were not present based on untreated file information."""

    # Fill out optionsDict with NA if true in untreated file
    untreatedSplit = untreatedName.split('_')
    if untreatedSplit[-2] == 'NA' and TYPE != 'LYSED':
        optionsDict['ANTIGENS HEALTHY'] = 'NA'

    if untreatedSplit[-3] == 'NA' and TYPE == 'LYSED':
        optionsDict['ANTIGENS HEALTHY'] = 'NA'

    return optionsDict

def save_data_subsetted(simsDF, saveLoc, xmlName, name, TYPE, states, subsetsRequested):
    """Save data subsetted into single pkl."""

    # Save file
    try:
        if states:
            with open(saveLoc + xmlName + '_' + name + 'STATES_' + TYPE + '.pkl', 'wb') as f:
                pickle.dump(simsDF, f)
        else:
            with open(saveLoc + xmlName + '_' + name + TYPE + '.pkl', 'wb') as f:
                pickle.dump(simsDF, f)

        scripts.subset.subset_utilities.print_save_message(subsetsRequested)

    except MemoryError:
        # Remove empty pickle file
        os.remove(saveLoc + xmlName + '_' + name + TYPE + '.pkl')
        print('MemoryError when attempting to save full pkl.')
        if TYPE == 'ANALYZED' and not states:
            print('Saving all non-list columns.')
            with open(saveLoc + xmlName + '_' + name + 'STATES_ANALYZED.pkl', 'wb') as f:
                simsDFstates = simsDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']
                pickle.dump(simsDFstates, f)
                # pickle.dump(simsDF[:, 0:-10], f)
            print('Saving all list columns in separate files.')
            simsDFcycle = simsDF[["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                                  "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "AVG CELL CYCLES CANCER",
                                  "AVG CELL CYCLES HEALTHY",
                                  "AVG CELL CYCLES T-CELL", "AVG CELL CYCLES CD4", "AVG CELL CYCLES CD8"]]
            simsDFvol = simsDF[["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                                "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "CELL VOLUMES CANCER",
                                "CELL VOLUMES HEALTHY",
                                "CELL VOLUMES T-CELL", "CELL VOLUMES CD4", "CELL VOLUMES CD8"]]
            with open(saveLoc + xmlName + '_' + name + 'CYCLES_ANALYZED.pkl', 'wb') as f:
                pickle.dump(simsDFcycle, f)
            with open(saveLoc + xmlName + '_' + name + 'VOLUMES_ANALYZED.pkl', 'wb') as f:
                pickle.dump(simsDFvol, f)

            scripts.subset.subset_utilities.print_save_message(subsetsRequested)

    return

def save_large_data_subsetted(simsDF, saveLoc, xmlName, name, TYPE, subsetsRequested):
    """Save large subsets into single pkl without list information."""

    print('Large number of files.')
    if TYPE == 'ANALYZED':
        print('Saving all non-list columns.')
        with open(saveLoc + xmlName + '_' + name + 'STATES_ANALYZED.pkl', 'wb') as f:
            pickle.dump(simsDF.iloc[:, 0:-10], f)
            scripts.subset.subset_utilities.print_save_message(subsetsRequested)

    return

def find_and_save_all_data(PKLFILES, files, simsDF, untreatedDF, states):
    """Find all data per file in list of files and add to datarame."""

    untreatedName = ''

    # Grab all files
    for file in PKLFILES:
        simDF = retrieve_file_data(file, states)
        if '0_NA_NA' in file:
            untreatedDF = untreatedDF.append(simDF, ignore_index=True)
            untreatedName = file.replace(files, '')
        else:
            simsDF = simsDF.append(simDF, ignore_index=True)

    return simsDF, untreatedDF, untreatedName

def collect_and_save_all_data(files, PKLFILES, xmlName, TYPE, saveLoc, states, subsetsRequested):
    """Collect all data in all files if no subsets selected."""

    simsDF, untreatedDF, TYPE = scripts.subset.subset_utilities.make_datatype_specific_df(TYPE)

    if (TYPE == 'ANALYZED' or TYPE == 'SHAREDLOCS') and states:
        simsDF, untreatedDF = drop_list_columns(simsDF, untreatedDF)

    # Grab all files
    simsDF, untreatedDF, untreatedName = find_and_save_all_data(PKLFILES, files, simsDF, untreatedDF, states)

    # Set up options dictionary for naming save file
    optionsDict = scripts.subset.subset_utilities.make_options_dict()

    # Fill out optionsDict with NA if true in untreated file
    optionsDict = update_options_dict_antigens_healthy_if_no_healthy_cells(optionsDict, untreatedName, TYPE)

    # Add untreated control to set
    simsDF = simsDF.append(untreatedDF, ignore_index=True)

    # Construct save file name
    name = scripts.subset.subset_utilities.construct_file_save_name(optionsDict)

    # Save file
    save_data_subsetted(simsDF, saveLoc, xmlName, name, TYPE, states, subsetsRequested)

    return

def find_and_save_files_within_specified_subset(PKLFILES, files, optionsDict, simsDF, untreatedDF, states):
    """Add files that belong in subset to dataframes."""

    print('Adding the following files to set:')

    fileCount = 0
    untreatedName = ''

    for file in PKLFILES:
        if '0_NA_NA' in file:
            print('\t' + file.replace(files, ''))
            simDF = retrieve_file_data(file, states)
            fileCount += 1
            untreatedDF = untreatedDF.append(simDF, ignore_index=True)
            untreatedName = file.replace(files, '')
        else:
            file_in_subset = check_if_file_in_subset(file, files, optionsDict)

            if file_in_subset:
                # Open file and save to set
                simDF = retrieve_file_data(file, states)
                fileCount += 1
                simsDF = simsDF.append(simDF, ignore_index=True)

    return simsDF, untreatedDF, untreatedName, fileCount

def collect_and_save_each_subset(subsetRequested, files, PKLFILES, xmlName, TYPE, saveLoc, states):
    """Collect all data in all files in given selected subset."""

    simsDF, untreatedDF, TYPE = scripts.subset.subset_utilities.make_datatype_specific_df(TYPE)

    if (TYPE == 'ANALYZED' or TYPE == 'SHAREDLOCS') and states:
        simsDF, untreatedDF = drop_list_columns(simsDF, untreatedDF)

    # Set up options dictionary for sorting files and naming save file
    optionsDict = scripts.subset.subset_utilities.make_options_dict()

    # Fill out options dictionary with set specifications
    for s in subsetRequested:
        optionsDict[s[0]] = str(s[1]).replace(':', '-')

    # Find only files with specificed set information
    simsDF, untreatedDF, untreatedName, fileCount = find_and_save_files_within_specified_subset(PKLFILES, files, optionsDict, simsDF, untreatedDF, states)

    # Fill out optionsDict with NA if true in untreated file
    optionsDict = update_options_dict_antigens_healthy_if_no_healthy_cells(optionsDict, untreatedName, TYPE)

    # Add untreated control to set
    simsDF = simsDF.append(untreatedDF, ignore_index=True)

    # Construct save file name
    name = scripts.subset.subset_utilities.construct_file_save_name(optionsDict)

    # Save set
    print('Saving set...')
    if (fileCount <= 6 and TYPE == 'ANALYZED') or (TYPE != 'ANALYZED'):
        save_data_subsetted(simsDF, saveLoc, xmlName, name, TYPE, states, subsetRequested)
    else:
        save_large_data_subsetted(simsDF, saveLoc, xmlName, name, TYPE, subsetRequested)

    return

def subset_data(files, xmlName, dataType, saveLoc, subsetsRequested='', states=False):
    """Collect desired subsets if any from list of all data and save in single file per subset."""

    # Parse requested subsets
    SUBSETS = parse_requested_subsets(subsetsRequested)

    # Get files
    PKLFILES = scripts.analyze.analyze_utilities.get_pkl_files(files)

    # Set type
    TYPE = dataType.upper()
    print(TYPE + ' SUBSET REQUESTED')

    # If no subsets specified, put all data in one file
    # Warning: Highly likely you will run into memory error
    if subsetsRequested == '':

        collect_and_save_all_data(files, PKLFILES, xmlName, TYPE, saveLoc, states, subsetsRequested)

    else:

        for subsetRequested in SUBSETS:

            collect_and_save_each_subset(subsetRequested, files, PKLFILES, xmlName, TYPE, saveLoc, states)

    return
