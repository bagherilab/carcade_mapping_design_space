import scripts.analyze.analyze_utilities
import scripts.plot.plot_cell_counts
import scripts.plot.plot_env
import scripts.plot.plot_heuristics
import scripts.plot.plot_kill_curves
import scripts.plot.plot_lysed
import scripts.plot.plot_spatial
import scripts.plot.plot_subcell_data
import scripts.plot.plot_utilities
import scripts.plot.plot_dish_tissue_compare
import pickle
import re
import pandas as pd
import matplotlib.pyplot as plt

def plot_dish_tissue_compare_data(files, color, saveLoc):
    """Call plotters that compare rank and score from dish and tissue simulations."""

    TREAT_RATIO_DICT_REVERSE = scripts.plot.plot_utilities.make_treat_ratio_reverse_key_dict()

    fileName = re.sub('.*RANK', 'RANK', files)
    print('\t' + fileName)

    with open(files, 'rb') as f:
        rankDF = pd.read_excel(f, header=0)

    rankDF['TREAT RATIO'] = rankDF['TREAT RATIO'].astype(str)
    rankDF['CAR AFFINITY'] = rankDF['CAR AFFINITY'].astype(str)
    for i in range(0, len(rankDF)):
        rankDF['TREAT RATIO'].iloc[i] = TREAT_RATIO_DICT_REVERSE[rankDF['TREAT RATIO'].iloc[i]]

    if color == 'X':
        for COLOR in ['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER']:
            scripts.plot.plot_dish_tissue_compare.plot_rank_ladder(rankDF, COLOR, 'VITRO', saveLoc)
            plt.close("all")

    return

def plot_heuristics_data(saveLoc):
    """Call plotters that plot binding heuristics."""

    scripts.plot.plot_heuristics.plot_binding_heuristic_CAR(saveLoc)
    scripts.plot.plot_heuristics.plot_binding_heuristic_self(saveLoc)

    return

def plot_kill_curve_exp_data(SAVELOC):
    """Call potters that plot experimental literature kill curve data."""

    scripts.plot.plot_kill_curves.plot_kill_curve_normalized_exp_separated(SAVELOC)

    return

def define_color_based_on_x_location(filesplit, X, COLOR, PARTIAL):
    """Define color based on X in file name, indicating all values of that feature present in file."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()

    if X == 1:
        if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X': COLOR = 'DOSE'
        if filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] == 'X': COLOR = 'TREAT RATIO'
        if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X': COLOR = 'CAR AFFINITY'
        if filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X': COLOR = 'ANTIGENS CANCER'
        if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'X': COLOR = 'ANTIGENS HEALTHY'

    if PARTIAL:
        COLOR = 'CAR AFFINITY'

    return COLOR

def plot_analyze_counts_data(simsDF, COLOR, ANALYSIS, filesplit, FILEID, fileid, SAVELOC):
    """Call plotters for plotting cell count dynamics for each cell population."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()

    print('\t\t' + 'Plotting count data over time')
    for p in range(0, len(POP_NAMES)):
        if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            scripts.plot.plot_cell_counts.plot_counts(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
            if ANALYSIS == 'ANALYZE':
                scripts.plot.plot_cell_counts.plot_counts_norm(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

            if 'VIVO' in FILEID and POP_NAMES[p] in ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE']:
                scripts.plot.plot_cell_counts.plot_counts_treat_norm(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

            if 'VIVO' in FILEID and POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                scripts.plot.plot_cell_counts.plot_counts_treat(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

            if 'LIVE' not in POP_NAMES[p]:
                scripts.plot.plot_cell_counts.plot_counts_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                if ANALYSIS == 'ANALYZE':
                    scripts.plot.plot_cell_counts.plot_counts_norm_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                if 'VIVO' in FILEID and POP_NAMES[p] in ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE']:
                    scripts.plot.plot_cell_counts.plot_counts_treat_norm_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                if 'VIVO' in FILEID and POP_NAMES[p] in ['CD4', 'CD8', 'T-CELL']:
                    scripts.plot.plot_cell_counts.plot_counts_treat_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

            plt.close("all")

    return

def plot_analyze_scatter_data(simsDF, COLOR, ANALYSIS, filesplit, fileid, SAVELOC):
    """Call plotters for analyzing living cancer vs healthy cell counts."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    # Plot scatter
    print('\t\t' + 'Plotting scatter data')
    if filesplit[FILEID_SPLIT_INDICES['POPS']] == 'CH':
        MARKER = 'o'
        if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
            scripts.plot.plot_cell_counts.plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, DISH_TIMES[-1])
            if ANALYSIS == 'ANALYZED':
                scripts.plot.plot_cell_counts.plot_CH_scatter_normalized(simsDF, COLOR, MARKER, fileid, SAVELOC, DISH_TIMES[-1])
        if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
            scripts.plot.plot_cell_counts.plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, TISSUE_TIMES[-1])
            if ANALYSIS == 'ANALYZED':
                scripts.plot.plot_cell_counts.plot_CH_scatter_normalized(simsDF, COLOR, MARKER, fileid, SAVELOC, TISSUE_TIMES[-1])
        plt.close("all")

    return

def plot_analyze_state_frac_data(simsDF, COLOR, fileid, SAVELOC):
    """Call plotters that plot cell state fractions over time."""

    # Plot type fracs
    print('\t\t' + 'Plotting state fraction data')
    scripts.plot.plot_subcell_data.plot_state_fracs(simsDF, COLOR, fileid, SAVELOC)
    scripts.plot.plot_subcell_data.plot_state_fracs(simsDF, COLOR, fileid, SAVELOC, True)
    scripts.plot.plot_subcell_data.plot_state_fracs_neutral(simsDF, COLOR, fileid, SAVELOC)
    if 'VIVO' in fileid:
        scripts.plot.plot_subcell_data.plot_state_fracs_treat(simsDF, COLOR, fileid, SAVELOC)
        scripts.plot.plot_subcell_data.plot_state_fracs_treat(simsDF, COLOR, fileid, SAVELOC, True)

    plt.close("all")

    return

def plot_analyze_volume_data(simsDF, COLOR, filesplit, FILEID, SAVELOC):
    """Call plotters that plot cell volume distributions."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    print('\t\t' + 'Plotting volume distribution data')
    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
        for time in range(1, len(DISH_TIMES)):
            scripts.plot.plot_subcell_data.plot_volumes(simsDF, COLOR, FILEID, SAVELOC, DISH_TIMES[time])
            scripts.plot.plot_subcell_data.plot_volumes_split(simsDF, COLOR, FILEID, SAVELOC, [4, 7])
        plt.close("all")
    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
        scripts.plot.plot_subcell_data.plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 1)
        scripts.plot.plot_subcell_data.plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 5)
        for time in TISSUE_TIMES:
            scripts.plot.plot_subcell_data.plot_volumes(simsDF, COLOR, FILEID, SAVELOC, time)
        plt.close("all")

    return

def plot_analyze_cycle_data(simsDF, COLOR, filesplit, FILEID, SAVELOC):
    """Call plotters that plot cell cycle distributions."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
        print('\t\t' + 'Plotting cycle distribution data')
        # Determine units: False = minutes, True = hours
        for unit in [False, True]:
            for time in range(1, len(DISH_TIMES)):
                scripts.plot.plot_subcell_data.plot_cycles(simsDF, COLOR, FILEID, SAVELOC, DISH_TIMES[time], unit)
            plt.close("all")
            scripts.plot.plot_subcell_data.plot_cycles_split(simsDF, COLOR, FILEID, SAVELOC, [4, 7], True)
        plt.close("all")

    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
        print('\t\t' + 'Plotting cycle distribution data')
        for unit in [False, True]:
            for time in TISSUE_TIMES:
                scripts.plot.plot_subcell_data.plot_cycles(simsDF, COLOR, FILEID, SAVELOC, time, unit)
            plt.close("all")
    return

def plot_kill_curve_sim_data(simsDF, filesplit, FILEID, fileid, SAVELOC):
    """Call plotters that plot simulated kill curve data."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    # Plot kill curve
    print('\t\t' + 'Plotting kill curve data')
    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
        for time in range(2, len(DISH_TIMES)):
            scripts.plot.plot_kill_curves.plot_kill_curve_normalized_sim(simsDF, FILEID, SAVELOC, 7)
    plt.close("all")

    return

def plot_env_conc_data(simsDF, COLOR, TIMES, FILEID, SAVELOC):
    """Call plotters that plot environment concentrations over time for each species."""

    MOL_NAMES = scripts.plot.plot_utilities.define_mol_names_list()

    print('\t\t' + 'Plotting environment concentration data')
    for m in range(0, len(MOL_NAMES)):
        if MOL_NAMES[m] == 'OXYGEN':
            continue
        else:
            scripts.plot.plot_env.plot_env_conc_times_bar(MOL_NAMES[m], simsDF, COLOR, FILEID, SAVELOC, TIMES)
            plt.close("all")

    return

def plot_env_parity_data(simsDF, COLOR, TIMES, FILEID, SAVELOC):
    """Call plotters that plot parity plots of environment species concentrations at the final time point in realistic vs ideal co-culture simulations."""

    MOL_NAMES = scripts.plot.plot_utilities.define_mol_names_list()

    print('\t\t' + 'Plotting environment parity data')
    for m in range(0, len(MOL_NAMES)):
        if MOL_NAMES[m] == 'OXYGEN':
            continue
        else:
            for XAXIS in ['IDEAL', 'REALISTIC']:
                if COLOR == 'X':
                    for color in ['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER']:
                        scripts.plot.plot_env.plot_evn_concs_ideal_realistic_parity(MOL_NAMES[m], simsDF, XAXIS, color, [TIMES[-1]], FILEID, SAVELOC)
                        plt.close("all")
                else:
                    scripts.plot.plot_env.plot_evn_concs_ideal_realistic_parity(MOL_NAMES[m], simsDF, XAXIS, COLOR, [TIMES[-1]], FILEID, SAVELOC)
                    plt.close("all")

    return

def plot_spatial_counts_data(simsDF, COLOR, TIMES, PARTIAL, filesplit, FILEID, SAVELOC):
    """Call plotters that plot cell spatial dynamics for each population"""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()

    print('\t\t' + 'Plotting count data across radius')
    for p in range(0, len(POP_NAMES)):
        for TIME in TIMES:
            if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
                continue
            else:
                scripts.plot.plot_spatial.plot_counts_radius(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC, TIME)
                scripts.plot.plot_spatial.plot_counts_radius(POP_NAMES[p] + ' NORMALIZED', simsDF, COLOR, FILEID, SAVELOC, TIME)
                plt.close("all")

            if 'LIVE' not in POP_NAMES[p]:
                scripts.plot.plot_spatial.plot_counts_radius_merge(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC, TIME)
                scripts.plot.plot_spatial.plot_counts_radius_merge(POP_NAMES[p] + ' NORMALIZED', simsDF, COLOR, FILEID, SAVELOC, TIME)

                plt.close("all")

    return

def plot_lysis_counts_data(simsDF, COLOR, filesplit, FILEID, SAVELOC):
    """Call plotters that plot cell lysis information for all lysed cell populations."""

    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    POP_LYSIS_NAMES = scripts.plot.plot_utilities.define_pop_lysis_names_list()

    # Plot counts data
    print('\t\t' + 'Plotting lysis count data over time')
    for p in range(0, len(POP_LYSIS_NAMES)):
        if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_LYSIS_NAMES[p]:
            continue
        else:
            scripts.plot.plot_lysed.plot_counts_lysed_time(POP_LYSIS_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
            scripts.plot.plot_lysed.plot_counts_lysed_time_exact(POP_LYSIS_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
        plt.close("all")

    if filesplit[FILEID_SPLIT_INDICES['POPS']] == 'CH':
        scripts.plot.plot_lysed.plot_counts_lysed_time_merge(POP_LYSIS_NAMES, simsDF, COLOR, FILEID, SAVELOC)
        scripts.plot.plot_lysed.plot_counts_lysed_time_exact_merge(POP_LYSIS_NAMES, simsDF, COLOR, FILEID, SAVELOC)
        plt.close("all")

    return

def plot_analyze_data(simsDF, COLOR, PARTIAL, ANALYSIS, FILEID, SAVELOC):
    """Plot and color data for analyze files based on file name."""

    filesplit = FILEID.split('_')
    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()

    # Count Xs in filesplit
    X = scripts.plot.plot_utilities.count_x_in_file_name(filesplit)

    if X == 1 or PARTIAL:
        if COLOR == 'X':
            COLOR = define_color_based_on_x_location(filesplit, X, COLOR, PARTIAL)

        fileid = FILEID.replace('_STATES','')

        # Plot counts data
        if 'CYCLES' not in FILEID and 'VOLUMES' not in FILEID:
            plot_analyze_counts_data(simsDF, COLOR, ANALYSIS, filesplit, FILEID, fileid, SAVELOC)
            plot_analyze_scatter_data(simsDF, COLOR, ANALYSIS, filesplit, fileid, SAVELOC)
            plot_analyze_state_frac_data(simsDF, COLOR, fileid, SAVELOC)

        # Plot distribution data
        if 'STATES' not in FILEID:
            # Plot volume distributions
            if 'CYCLES' not in FILEID:
                plot_analyze_volume_data(simsDF, COLOR, filesplit, FILEID, SAVELOC)

            # Plot cycle distributions
            if 'VOLUMES' not in FILEID:
                plot_analyze_cycle_data(simsDF, COLOR, filesplit, FILEID, SAVELOC)

    if 'VOLUMES' not in FILEID and 'CYCLES' not in FILEID and not PARTIAL:
        fileid = FILEID.replace('_STATES', '')
        if filesplit[FILEID_SPLIT_INDICES['DOSE']] != 'X' or filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] != 'X':
            if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X' and filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X':
                if COLOR == 'X': COLOR = 'CAR AFFINITY'

                # Plot kill curve
                plot_kill_curve_sim_data(simsDF, filesplit, FILEID, fileid, SAVELOC)

    return

def plot_env_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC):
    """Plot and color data for environment files based on file name."""

    filesplit = FILEID.split('_')
    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    # Count Xs in filesplit
    X = scripts.plot.plot_utilities.count_x_in_file_name(filesplit)

    if filesplit[FILEID_SPLIT_INDICES['EXP']] == 'VITRO':
        TIMES = [i for i in DISH_TIMES if i != 0]
    else:
        TIMES = TISSUE_TIMES

    if X == 1 or PARTIAL:
        if COLOR == 'X':
            COLOR = define_color_based_on_x_location(filesplit, X, COLOR, PARTIAL)

        plot_env_conc_data(simsDF, COLOR, TIMES, FILEID, SAVELOC)

    if X == 5 and 'VITRO' in FILEID:
        plot_env_parity_data(simsDF, COLOR, TIMES, FILEID, SAVELOC)

    return

def plot_spatial_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC):
    """Plot and color data for spatial files based on file name."""

    filesplit = FILEID.split('_')
    FILEID_SPLIT_INDICES = scripts.plot.plot_utilities.define_fileid_split_indices_dict()
    DISH_TIMES = scripts.plot.plot_utilities.define_dish_times_list()
    TISSUE_TIMES = scripts.plot.plot_utilities.define_tissue_times_list()

    if filesplit[FILEID_SPLIT_INDICES['EXP']] == 'VITRO':
        TIMES = DISH_TIMES
    else:
        TIMES = TISSUE_TIMES

    # Count Xs in filesplit
    X = scripts.plot.plot_utilities.count_x_in_file_name(filesplit)

    if X == 1 or PARTIAL:
        if COLOR == 'X':
            COLOR = define_color_based_on_x_location(filesplit, X, COLOR, PARTIAL)

        # Plot counts data
        plot_spatial_counts_data(simsDF, COLOR, TIMES, PARTIAL, filesplit, FILEID, SAVELOC)

    return

def plot_lysed_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC):
    """Plot and color data for lysed files based on file name."""

    filesplit = FILEID.split('_')

    # Count Xs in filesplit
    X = scripts.plot.plot_utilities.count_x_in_file_name(filesplit)

    if X == 1 or PARTIAL:
        if COLOR == 'X':
            COLOR = define_color_based_on_x_location(filesplit, X, COLOR, PARTIAL)

    plot_lysis_counts_data(simsDF, COLOR, filesplit, FILEID, SAVELOC)

    return

def determine_data_file_type(fileName):
    """Determine data file type (and thus plots to make) based on file name."""

    if 'ANALYZED' in fileName:
        FILEID = fileName.replace('_ANALYZED.pkl', '')
        analysis = 'ANALYZED'

    if 'ENVIRONMENT' in fileName:
        FILEID = fileName.replace('_ENVIRONMENT.pkl', '')
        analysis = 'ENVIRONMENT'

    if 'SPATIAL' in fileName:
        FILEID = fileName.replace('_SPATIAL.pkl', '')
        analysis = 'SPATIAL'

    if 'LYSED' in fileName:
        FILEID = fileName.replace('_LYSED.pkl', '')
        analysis = 'LYSED'

    if 'SHAREDLOCS' in fileName:
        FILEID = fileName.replace('_SHAREDLOCS.pkl', '')
        analysis = 'SHAREDLOCS'

    return FILEID, analysis

def plot_data_based_on_analysis_type(simsDF, ANALYSIS, COLOR, PARTIAL, FILEID, SAVELOC):
    """Call appropriate data plotter sequence based on file type."""

    if ANALYSIS == 'ANALYZED' or ANALYSIS == 'SHAREDLOCS':

        plot_analyze_data(simsDF, COLOR, PARTIAL, ANALYSIS, FILEID, SAVELOC)

    if ANALYSIS == 'ENVIRONMENT':

        plot_env_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC)

    if ANALYSIS == 'SPATIAL':

        plot_spatial_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC)

    if ANALYSIS == 'LYSED':

        plot_lysed_data(simsDF, COLOR, PARTIAL, FILEID, SAVELOC)

    return

def plot_data(files, color, saveLoc='', partial=False):
    """Iterate through all files and plot appropriate file types.

    plot_data takes a directory of (or a single) .pkl simulation files that result from analyze_cells, analyze_env, analyze_spatial, or analyze_lysis and
    plots the data features in the dataframe over time.

    Usage:
        plot_data(files, color, saveLoc='', partial=False)

        files
            Path to .pkl or directory.
        color
            Feature by which to color data by. If X given, will color by axis along which data varies.
            If two featuers vary, default will be used.
        [partial]
            Flag indicating only partial dataset present instead of full combinatoral set.
        [saveLoc]
            Location of where to save file, default will save here.
    """

    print("Making figures for the following files:")

    # Get files
    PKLFILES = scripts.analyze.analyze_utilities.get_pkl_files(files)

    for file in PKLFILES:

        fileName = scripts.plot.plot_utilities.get_file_name(file)

        print('\t' + fileName)

        FILEID, analysis = determine_data_file_type(fileName)

        with open(file, 'rb') as f:
            simsDF = pickle.load(f)

        plot_data_based_on_analysis_type(simsDF, analysis, color, partial, FILEID, saveLoc)

    print("Finished making plots for all files.")

    return