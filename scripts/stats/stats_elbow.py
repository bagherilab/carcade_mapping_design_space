import scripts.plot.plot_utilities
import scripts.stats.stats_utilities
import matplotlib.pyplot as plt

def plot_output_elbow_sort_output_and_color(simsDF, OUTPUT, COLOR, FILEID, NORM, SCORE, SAVELOC):
    """Plot elbow plots and sort by output and color by feature selected."""

    COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = [OUTPUT] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i][OUTPUT], color=colorDict[key])

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(OUTPUT + " VALUE", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(OUTPUT + " VALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_sort_all(simsDF, OUTPUT, COLOR, FILEID, NORM, SCORE, SAVELOC):
    """Plot elbow plots and sort by  by features in subsequent order."""

    COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    sortby = [OUTPUT] + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i][OUTPUT], color=colorDict[key])

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(OUTPUT + "\nVALUE", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(OUTPUT + "\nVALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = '_' + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_ALLSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

def plot_output_elbow_CH_sort_cancer_and_color(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):
    """Plot elbow plots of normalized cancer and healthy values on same plot and sort by normalized cancer value and color by feature selected."""

    COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = ['Y_NORM_CANCER_LIVE'] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = '_' + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_CH_CANCER_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_CH_sort_healthy_and_color(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):
    """Plot elbow plots of normalized cancer and healthy values on same plot and sort by normalized healthy value and color by feature selected."""

    COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = ['Y_NORM_HEALTHY_LIVE'] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_CH_HEALTHY_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_CH_sort_all(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):
    """Plot elbow plots of normalized cancer and healthy values on same plot and sort by all outputs sequentially and color by feature selected."""

    COLOR_DICT_HEATMAP = scripts.stats.stats_utilities.make_features_heatmap_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    sortby = ['Y_NORM_HEALTHY_LIVE', 'Y_NORM_CANCER_LIVE'] + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_ALLSORT_CH_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')