import scripts.plot.plot_utilities
import scripts.stats.stats_utilities
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors as mcolors

def check_for_extended_dose_data(simsDF, X, Y):
    countExtended = 0

    for i in range(0, len(simsDF)):
        value = simsDF.iloc[i]['DOSE']
        if value == 5000 or value == 10000:
            if X == 'DOSE':
                countExtended += 1
                break
            if Y == 'DOSE':
                countExtended += 1
                break

    if countExtended == 0:

        if X == 'DOSE':
            X_AXIS = [250, 500, 1000]

            return X_AXIS

        if Y == 'DOSE':
            Y_AXIS = [250, 500, 1000]

            return Y_AXIS

def check_for_extended_treat_ratio_data(simsDF, X, Y):

    countExtended = 0

    for i in range(0,len(simsDF)):
        value = simsDF.iloc[i]['TREAT_RATIO']
        if value == '10:90' or value == '90:10':
            if X == 'TREAT RATIO':
                countExtended += 1
                break
            if Y == 'TREAT RATIO':
                countExtended += 1
                break

    if countExtended == 0:

        if X == 'TREAT RATIO':
            X_AXIS = ["0:100", "25:75", "50:50", "75:25", "100:0"]

            return X_AXIS

        if Y == 'TREAT RATIO':
            Y_AXIS = ["0:100", "25:75", "50:50", "75:25", "100:0"]

            return Y_AXIS

def plot_two_axis_heatmaps_diverging_all(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):
    """Plot two axes heatmaps with diverging color map where all possible values of all features not shown are included."""

    AXES_SETS = scripts.stats.stats_utilities.define_axes_sets_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if X == 'TREAT RATIO' or Y == 'TREAT RATIO':
        if X == 'TREAT RATIO':
            X_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)

    if X == 'DOSE' or Y == 'DOSE':
        if X == 'DOSE':
            X_AXIS = check_for_extended_dose_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_dose_data(simsDF, X, Y)

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    response = 'CANCER'
    cmap = 'coolwarm'
    center = 1.0
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
        vmin = 0
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        cmap = 'coolwarm_r'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
        cmap = 'coolwarm_r'
    if COLOR == 'SCORE':
        response = 'SCORE'
        cmap = 'coolwarm_r'
        vmax = 2
        center = 0

    heatmapDict = {}

    for y in Y_AXIS:
        heatmapDict[y] = {}
        for x in X_AXIS:
            heatmapDict[y][x] = []

    for y in Y_AXIS:
        for x in X_AXIS:
            if Y == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
            elif X == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
            else:
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

            for i in range(0,len(subsetDF)):
                heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

    heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

    for y in Y_AXIS:
        for x in X_AXIS:
            heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

    figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figHeatmap.add_subplot(1, 1, 1)
    ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                     xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                     cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
    ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(COLOR + "\nRESPONSE - ALL", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    if ANNOTATE:
        annotate = '_ANNOTATED'
    else:
        annotate = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_DIVERGING_ALL_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_diverging_constant(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):
    """Plot two axes heatmaps with diverging color map where features not shown are held at intermediate values."""

    AXES_SETS = scripts.stats.stats_utilities.define_axes_sets_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    CONSTANTS = scripts.stats.stats_utilities.define_constants_dict()

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if X == 'TREAT RATIO' or Y == 'TREAT RATIO':
        if X == 'TREAT RATIO':
            X_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)

    if X == 'DOSE' or Y == 'DOSE':
        if X == 'DOSE':
            X_AXIS = check_for_extended_dose_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_dose_data(simsDF, X, Y)

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    for key in AXES_SETS:
        key_ = key.replace(" ",'_')
        if key != X and key != Y and key != 'ANTIGENS HEALTHY':
            simsDF = simsDF[simsDF[key_] == CONSTANTS[key]]

    DFlist = []

    if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
        for i in range(0, len(CONSTANTS['ANTIGENS HEALTHY'])):
            DFlist.append(simsDF[simsDF['ANTIGENS_HEALTHY'] == CONSTANTS['ANTIGENS HEALTHY'][i]])
    else:
        DFlist.append(simsDF)

    response = 'CANCER'
    cmap = 'coolwarm'
    center = 1.0
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
        vmin = 0
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        cmap = 'coolwarm_r'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        cmap = 'coolwarm_r'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        cmap = 'coolwarm_r'
        vmax = 2
        center = 0

    for d in range(0,len(DFlist)):

        DF = DFlist[d]
        heatmapDict = {}

        for y in Y_AXIS:
            heatmapDict[y] = {}
            for x in X_AXIS:
                heatmapDict[y][x] = []

        for y in Y_AXIS:
            for x in X_AXIS:
                if Y == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
                elif X == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
                else:
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

                for i in range(0, len(subsetDF)):
                    heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

        heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

        for y in Y_AXIS:
            for x in X_AXIS:
                heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

        figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
        ax = figHeatmap.add_subplot(1, 1, 1)
        ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                         xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                         cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
        ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_title(COLOR + "\nRESPONSE - CONSTANT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        add_HA = ''
        if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
            add_HA = '_HA_' + str(CONSTANTS['ANTIGENS HEALTHY'][d])
        ax.set_title(COLOR + " RESPONSE - CONSTANT" + add_HA.replace('_',''), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

        if '_CH_' in FILEID:
            add_score = "_" + SCORE
        else:
            add_score = ''

        if ANNOTATE:
            annotate = '_ANNOTATED'
        else:
            annotate = ''

        plt.xticks(fontsize=TICKSIZE)
        plt.yticks(fontsize=TICKSIZE)

        plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_CONSTANT_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + add_HA + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_grayscale_all(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):
    """Plot two axes heatmaps with diverging color map where all possible values of all features not shown are included."""

    AXES_SETS = scripts.stats.stats_utilities.define_axes_sets_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if X == 'TREAT RATIO' or Y == 'TREAT RATIO':
        if X == 'TREAT RATIO':
            X_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)

    if X == 'DOSE' or Y == 'DOSE':
        if X == 'DOSE':
            X_AXIS = check_for_extended_dose_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_dose_data(simsDF, X, Y)

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    response = 'CANCER'
    cmap = 'Greys'
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        vmax = 2

    heatmapDict = {}

    for y in Y_AXIS:
        heatmapDict[y] = {}
        for x in X_AXIS:
            heatmapDict[y][x] = []

    for y in Y_AXIS:
        for x in X_AXIS:
            if Y == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
            elif X == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
            else:
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

            for i in range(0,len(subsetDF)):
                heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

    heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

    for y in Y_AXIS:
        for x in X_AXIS:
            heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

    figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figHeatmap.add_subplot(1, 1, 1)
    ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                     xticklabels=X_AXIS, yticklabels=Y_AXIS, cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
    ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(COLOR + "\nRESPONSE - ALL", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    if ANNOTATE:
        annotate = '_ANNOTATED'
    else:
        annotate = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_GRAYSCALE_ALL_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_grayscale_constant(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):
    """Plot two axes heatmaps with grayscale color map where features not shown are held at intermediate values."""

    AXES_SETS = scripts.stats.stats_utilities.define_axes_sets_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    CONSTANTS = scripts.stats.stats_utilities.define_constants_dict()

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if X == 'TREAT RATIO' or Y == 'TREAT RATIO':
        if X == 'TREAT RATIO':
            X_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_treat_ratio_data(simsDF, X, Y)

    if X == 'DOSE' or Y == 'DOSE':
        if X == 'DOSE':
            X_AXIS = check_for_extended_dose_data(simsDF, X, Y)
        else:
            Y_AXIS = check_for_extended_dose_data(simsDF, X, Y)

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    for key in AXES_SETS:
        key_ = key.replace(" ",'_')
        if key != X and key != Y and key != 'ANTIGENS HEALTHY':
            simsDF = simsDF[simsDF[key_] == CONSTANTS[key]]

    DFlist = []

    if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
        for i in range(0, len(CONSTANTS['ANTIGENS HEALTHY'])):
            DFlist.append(simsDF[simsDF['ANTIGENS_HEALTHY'] == CONSTANTS['ANTIGENS HEALTHY'][i]])
    else:
        DFlist.append(simsDF)

    response = 'CANCER'
    cmap = 'Greys'
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        vmax = 2

    for d in range(0,len(DFlist)):

        DF = DFlist[d]
        heatmapDict = {}

        for y in Y_AXIS:
            heatmapDict[y] = {}
            for x in X_AXIS:
                heatmapDict[y][x] = []

        for y in Y_AXIS:
            for x in X_AXIS:
                if Y == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
                elif X == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
                else:
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

                for i in range(0, len(subsetDF)):
                    heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

        heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

        for y in Y_AXIS:
            for x in X_AXIS:
                heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

        figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
        ax = figHeatmap.add_subplot(1, 1, 1)
        ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                         xticklabels=X_AXIS, yticklabels=Y_AXIS, cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
        ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_title(COLOR + "\nRESPONSE - CONSTANT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        add_HA = ''
        if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
            add_HA = '_HA_' + str(CONSTANTS['ANTIGENS HEALTHY'][d])
        ax.set_title(COLOR + " RESPONSE - CONSTANT" + add_HA.replace('_',''), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

        if '_CH_' in FILEID:
            add_score = "_" + SCORE
        else:
            add_score = ''

        if ANNOTATE:
            annotate = '_ANNOTATED'
        else:
            annotate = ''

        plt.xticks(fontsize=TICKSIZE)
        plt.yticks(fontsize=TICKSIZE)

        plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_CONSTANT_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + add_HA + annotate + '.svg', bbox_inches='tight')

    return

def plot_output_heatmap_diverging(simsDF, OUTPUT, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with chosen output shown with diverging color map."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1.0)
            plt.register_cmap(name=cm_name, data=cm)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    ascending = []
    for c in columns:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_CANCER_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=columns, ascending=ascending)

    columns.append(OUTPUT)
    if 'CANCER' in OUTPUT:
        cmap_response = 'coolwarm'
        center = 1.0
        vmax = 3
        if '_CH_' in FILEID:
            vmax = 5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 0.35
    elif 'HEALTHY' in OUTPUT:
        cmap_response = 'coolwarm_r'
        center = 1.0
        vmax = 1.5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 1
    elif 'TCELL' in OUTPUT:
        cmap_response = 'coolwarm_r'
        center = 1.0
        vmax = 300
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 100
    else:
        cmap_response = 'coolwarm_r'
        center = 0
        vmax = 2
    cmap.append(cmap_response)
    vmaxes.append(vmax)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0})

    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=True, center=center, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_DIVERGING_" + OUTPUT + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale(simsDF, OUTPUT, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with chosen output shown with grayscale color map."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1.0)
            plt.register_cmap(name=cm_name, data=cm)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    ascending = []
    for c in columns:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_CANCER_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=columns, ascending=ascending)

    columns.append(OUTPUT)
    cmap_response = 'Greys'
    if 'CANCER' in OUTPUT:
        vmax = 3
        if '_CH_' in FILEID:
            vmax = 5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 0.35
    elif 'HEALTHY' in OUTPUT:
        vmax = 1.5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 1
    elif 'TCELL' in OUTPUT:
        vmax = 300
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 100
    else:
        vmax = 2
    cmap.append(cmap_response)
    vmaxes.append(vmax)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0})

    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=True, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_GRAYSCALE_" + OUTPUT + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_diverging_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with all possible outputs shown with diverging color map and optionally sorted by score."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    if SORT_SCORES: sortby = ['SCORE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=m, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_SCORES' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with all possible outputs shown with grayscale color map and optionally sorted by score."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    if SORT_SCORES: sortby = ['SCORE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    cmap_response = 'Greys'
    for output in OUTPUTS:
        if 'CANCER' in output:
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_SCORES' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDF, OUTPUTS, SORT, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with all possible outputs shown with diverging color map and sorted by chosen output."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    sortby = [SORT] + columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=m, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_'+SORT
    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDF, OUTPUTS, SORT, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with all possible outputs shown with grayscale color map and sorted by chosen output."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)

    sortby = [SORT] + columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    cmap_response = 'Greys'
    for output in OUTPUTS:
        if 'CANCER' in output:
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_'+SORT
    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_with_subplots_line_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    """Plot heatmap of all feature values with all possible outputs shown with line plots and sorted by output based on simulation type."""

    if max(simsDF['DOSE']) > 1000:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict_extended_dose()
    else:
        COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []
    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = scripts.stats.stats_utilities.NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = scripts.stats.stats_utilities.maximum_absolute_scaling(simsDF)
    # print(simsDF['ANTIGENS HEALTHY'])

    if SORT_SCORES:
        if '_CH_' in FILEID:
            sortby = ['SCORE'] + columns
        else:
            sortby = ['Y_NORM_CANCER_LIVE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(True)
        else:
            ascending.append(False)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(len(columns),1, gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            x = np.array([simsDF[columns[i]].values])
            sns.heatmap(x, ax=a, xticklabels=[], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        elif 'SCORE' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            y = [0 for j in x]
            sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i-1], color='black')
            sns.lineplot(x, y, ax=ax[i-1], color='black')
            ax[i-1].set_xlim([0,len(simsDF)-1])
            ax[i-1].set_ylim([min(simsDF[columns[i]]),max(simsDF['SCORE'])])
            ax[i-1].set_yticks([min(simsDF[columns[i]]),max(simsDF[columns[i]])])
            ax[i-1].set_yticklabels([min(simsDF[columns[i]]),max(simsDF[columns[i]])])
            a.set_visible(False)
            ax[i-1].set_ylabel('SCORE', rotation='horizontal', labelpad=15)
        elif 'CANCER' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            sns.lineplot(x, simsDF[columns[i]].values, ax=a, color='darkgray')
            a.set_xlim([0, len(simsDF)-1])
            a.set_ylim([0, max(simsDF[columns[i]])])
            a.set_yticks([0,max(simsDF[columns[i]])])
            a.set_yticklabels([0,max(simsDF[columns[i]])])
            a.set_ylabel('NORMALIZED CANCER AND HEALTHY VALUE',rotation='horizontal', labelpad=125)
            if '_C_' in FILEID:
                y = [1 for i in x]
                sns.lineplot(x, y, ax=a, color='lightgray')
            if '_CH_' in FILEID:
                y = [1 for i in x]
                sns.lineplot(x, y, ax=a, color='black', zorder=0)
        elif 'HEALTHY' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i-1], color='lightgray')
            if max(simsDF[columns[i]]) > max(simsDF[columns[i-1]]):
                ax[i-1].set_ylim([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticks([0,max(simsDF[columns[i]])])
                ax[i-1].set_yticklabels([0, max(simsDF[columns[i]])])
        else:
            x = [index for index in range(0, len(simsDF))]
            if '_CH_' in FILEID:
                sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i - 1], color='darkgray')
                ax[i-1].set_xlim([0, len(simsDF)-1])
                ax[i-1].set_ylim([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticks([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticklabels([0, max(simsDF[columns[i]])])
                ax[i-1].set_ylabel(columns[i], rotation='horizontal', labelpad=75)
            else:
                sns.lineplot(x, simsDF[columns[i]].values, ax=a, color='darkgray')
                a.set_xlim([0, len(simsDF)-1])
                a.set_ylim([0, max(simsDF[columns[i]])])
                a.set_yticks([0, max(simsDF[columns[i]])])
                a.set_yticklabels([0, max(simsDF[columns[i]])])
                a.set_ylabel(columns[i], rotation='horizontal', labelpad=75)
        if i < 5 and '_CH_' in FILEID:
            a.set_yticklabels([s], rotation=0)
        elif i < 4 and '_C_' in FILEID:
            a.set_yticklabels([s], rotation=0)
        a.xaxis.set_ticks([])


    if '_CH_' in FILEID:
        add_sort = '_SORTED_SCORES' if SORT_SCORES else ''
    else:
        add_sort = '_SORTED_CANCER' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_ALLOUTPUTS_SUBPLOTS" + add_sort + '.pdf', bbox_inches='tight')

    return