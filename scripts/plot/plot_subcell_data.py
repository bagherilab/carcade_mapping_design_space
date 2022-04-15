import scripts.analyze.analyze_utilities
import scripts.plot.plot_utilities
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_state_fracs(simsDF, COLOR, FILEID, SAVELOC):
    """Plot fractions of cells in each state over time for all cell populations and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figStateFracs, axs = plt.subplots(12, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):

        if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = 0
        elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = "CONTROL"
        else:
            color = simsDF.iloc[i][COLOR]

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            start = 0
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            start = 2

        # Plot CANCER cell state fracs
        axs[0, 0].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 0].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 0].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 0].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 0].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 0].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 0].set_visible(False)
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 1].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 1].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 1].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 1].set_visible(False)
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 2].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 2].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 2].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 2].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 2].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 2].set_visible(False)
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 3].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 3].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 3].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 3].set_visible(False)
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 4].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 4].set_visible(False)
        axs[3, 4].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 4].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 4].set_visible(False)
        axs[6, 4].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 4].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 4].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 4].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 4].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 5].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 5].set_visible(False)
        axs[3, 5].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 5].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 5].set_visible(False)
        axs[6, 5].set_visible(False)
        axs[7, 5].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 5].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 5].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 5].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 6].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 6].set_visible(False)
        axs[3, 6].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 6].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 6].set_visible(False)
        axs[6, 6].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 6].set_visible(False)
        axs[8, 6].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 6].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 6].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 6].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])

    for a in range(0, 12):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(plot_time)])
            axs[a, b].set_ylim([0, 1])

    YLABELS = ['MIGRAGORY', 'PROLIFERATIVE', 'QUIESCENT', 'SENESCENT', 'APOPTOTIC',
               'NECROTIC', 'CYTOTOXIC', 'STIMULATORY', 'EXHAUSTED', 'ANERGIC', 'STARVED', 'PAUSED']
    TITLES = ['CANCER CELLS', 'ALIVE CANCER CELLS', 'HEALTHY CELLS', 'ALIVE HEALTHY CELLS',
              'CAR T-CELLS', 'CD4 CAR T-CELLS', 'CD8 CAR T-CELLS']

    listYAxes = [0, 7, 14, 21, 28, 30, 32, 35, 37, 46, 48, 53, 60, 67, 74, 81]
    listXAxes = [11, 12, 13, 22, 24, 32, 33, 34, 35, 37, 48, 81, 82, 83]
    listYlabels = [0, 7, 14, 21, 28, 35, 46, 53, 60, 67, 74, 81]
    j = 0
    k = 0
    for i, ax in enumerate(axs.flat):
        if i in listYAxes and i in listYlabels:
            ax.set_ylabel(YLABELS[j], fontname='Arial', fontweight='bold', fontsize=7, labelpad=10)
            ax.get_yaxis().set_visible(True)
            j = j + 1
        elif i in listYAxes and i not in listYlabels:
            ax.get_yaxis().set_visible(True)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 7:
            ax.set_title(TITLES[k], fontname='Arial', fontweight='bold', fontsize=7, pad=10)
            k = k + 1
        if i in listXAxes:
            ax.get_xaxis().set_visible(True)
        else:
            ax.get_xaxis().set_visible(False)

    figStateFracs.text(0.5, 0.08, "TIME", ha='center', fontname='Arial', fontweight='bold', fontsize=13)
    figStateFracs.text(0.01, 0.5, "STATE (FRACTION)", va='center', rotation='vertical', fontname='Arial',
                     fontweight='bold', fontsize=13)

    for color in colorDict:
        axs[0,0].plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        axs[0,0].plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    axs[0, 0].legend(bbox_to_anchor=(9.5, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_STATESFRAC.svg', bbox_inches='tight')

    return

def plot_state_fracs_treat(simsDF, COLOR, FILEID, SAVELOC):
    """Plot fractions of cells in each state over treatment time for all cell populations and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figStateFracs, axs = plt.subplots(12, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):

        if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = 0
        elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = "CONTROL"
        else:
            color = simsDF.iloc[i][COLOR]

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            start = 0
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][44:]]
            start = 44

        # Plot CANCER cell state fracs
        axs[0, 0].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 0].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 0].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 0].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 0].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 0].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 0].set_visible(False)
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 1].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 1].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 1].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 1].set_visible(False)
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 2].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 2].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 2].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 2].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 2].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 2].set_visible(False)
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 3].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 3].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 3].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 3].set_visible(False)
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 4].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 4].set_visible(False)
        axs[3, 4].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 4].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 4].set_visible(False)
        axs[6, 4].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 4].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 4].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 4].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 4].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 5].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 5].set_visible(False)
        axs[3, 5].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 5].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 5].set_visible(False)
        axs[6, 5].set_visible(False)
        axs[7, 5].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 5].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 5].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 5].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 6].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 6].set_visible(False)
        axs[3, 6].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 6].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 6].set_visible(False)
        axs[6, 6].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 6].set_visible(False)
        axs[8, 6].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 6].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 6].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 6].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])

    for a in range(0, 12):
        for b in range(0, 7):
            axs[a, b].set_xlim([21, max(plot_time)])
            axs[a, b].set_ylim([0, 1])

    YLABELS = ['MIGRAGORY', 'PROLIFERATIVE', 'QUIESCENT', 'SENESCENT', 'APOPTOTIC',
               'NECROTIC', 'CYTOTOXIC', 'STIMULATORY', 'EXHAUSTED', 'ANERGIC', 'STARVED', 'PAUSED']
    TITLES = ['CANCER CELLS', 'ALIVE CANCER CELLS', 'HEALTHY CELLS', 'ALIVE HEALTHY CELLS',
              'CAR T-CELLS', 'CD4 CAR T-CELLS', 'CD8 CAR T-CELLS']

    listYAxes = [0, 7, 14, 21, 28, 30, 32, 35, 37, 46, 48, 53, 60, 67, 74, 81]
    listXAxes = [11, 12, 13, 22, 24, 32, 33, 34, 35, 37, 48, 81, 82, 83]
    listYlabels = [0, 7, 14, 21, 28, 35, 46, 53, 60, 67, 74, 81]
    j = 0
    k = 0
    for i, ax in enumerate(axs.flat):
        if i in listYAxes and i in listYlabels:
            ax.set_ylabel(YLABELS[j], fontname='Arial', fontweight='bold', fontsize=7, labelpad=10)
            ax.get_yaxis().set_visible(True)
            j = j + 1
        elif i in listYAxes and i not in listYlabels:
            ax.get_yaxis().set_visible(True)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 7:
            ax.set_title(TITLES[k], fontname='Arial', fontweight='bold', fontsize=7, pad=10)
            k = k + 1
        if i in listXAxes:
            ax.get_xaxis().set_visible(True)
        else:
            ax.get_xaxis().set_visible(False)

    figStateFracs.text(0.5, 0.08, "TIME", ha='center', fontname='Arial', fontweight='bold', fontsize=13)
    figStateFracs.text(0.01, 0.5, "STATE (FRACTION)", va='center', rotation='vertical', fontname='Arial',
                     fontweight='bold', fontsize=13)

    for color in colorDict:
        axs[0,0].plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        axs[0,0].plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    axs[0, 0].legend(bbox_to_anchor=(9.5, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_STATESFRAC_TREAT.svg', bbox_inches='tight')

    return

def plot_state_fracs_neutral(simsDF, COLOR, FILEID, SAVELOC):
    """Plot fractions of cells in each state (including neutral state) over time for all cell populations  and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figStateFracs, axs = plt.subplots(13, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):

        if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = 0
        elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = "CONTROL"
        else:
            color = simsDF.iloc[i][COLOR]

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            start = 0
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            start = 2


        # Plot CANCER cell state fracs
        axs[0, 0].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 0].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 0].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 0].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 0].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 0].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 0].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)
        axs[12, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 1].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 1].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 1].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 1].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)
        axs[12, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 2].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 2].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 2].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 2].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 2].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 2].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)
        axs[12, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 3].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 3].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 3].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 3].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)
        axs[12, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 4].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 4].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 4].set_visible(False)
        axs[4, 4].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 4].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 4].set_visible(False)
        axs[7, 4].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 4].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 4].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 4].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 5].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 5].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 5].set_visible(False)
        axs[4, 5].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 5].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 5].set_visible(False)
        axs[7, 5].set_visible(False)
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 5].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 5].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 5].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 6].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 6].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 6].set_visible(False)
        axs[4, 6].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 6].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 6].set_visible(False)
        axs[7, 6].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 6].set_visible(False)
        axs[9, 6].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 6].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 6].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 6].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])

    for a in range(0, 13):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(plot_time)])
            axs[a, b].set_ylim([0, 1])

    YLABELS = ['NEUTRAL', 'MIGRAGORY', 'PROLIFERATIVE', 'QUIESCENT', 'SENESCENT', 'APOPTOTIC',
               'NECROTIC', 'CYTOTOXIC', 'STIMULATORY', 'EXHAUSTED', 'ANERGIC', 'STARVED', 'PAUSED']
    TITLES = ['CANCER CELLS', 'ALIVE CANCER CELLS', 'HEALTHY CELLS', 'ALIVE HEALTHY CELLS',
              'CAR T-CELLS', 'CD4 CAR T-CELLS', 'CD8 CAR T-CELLS']

    listYAxes = [0, 7, 14, 21, 28, 35, 37, 39, 42, 44, 53, 55, 60, 67, 74, 81, 88]
    listXAxes = [18, 19, 20, 29, 31, 39, 40, 41, 42, 44, 55, 88, 89, 90]
    listYlabels = [0, 7, 14, 21, 28, 35, 42, 53, 60, 67, 74, 81, 88]
    j = 0
    k = 0
    for i, ax in enumerate(axs.flat):
        if i in listYAxes and i in listYlabels:
            ax.set_ylabel(YLABELS[j], fontname='Arial', fontweight='bold', fontsize=7, labelpad=10)
            ax.get_yaxis().set_visible(True)
            j = j + 1
        elif i in listYAxes and i not in listYlabels:
            ax.get_yaxis().set_visible(True)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 7:
            ax.set_title(TITLES[k], fontname='Arial', fontweight='bold', fontsize=7, pad=10)
            k = k + 1
        if i in listXAxes:
            ax.get_xaxis().set_visible(True)
        else:
            ax.get_xaxis().set_visible(False)

    figStateFracs.text(0.5, 0.08, "TIME", ha='center', fontname='Arial', fontweight='bold', fontsize=13)
    figStateFracs.text(0.01, 0.5, "STATE (FRACTION)", va='center', rotation='vertical', fontname='Arial',
                     fontweight='bold', fontsize=13)

    for color in colorDict:
        axs[0,0].plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        axs[0,0].plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    axs[0, 0].legend(bbox_to_anchor=(9.5, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_STATESFRACNEUTR.svg', bbox_inches='tight')

    return

def plot_volumes(simsDF, COLOR, FILEID, SAVELOC, TIME):
    """Plot cell volume distributions at a given time point for each cell population and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:
                key = 'CELL VOLUMES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []

                for i in range(0, len(simsDF)):
                    y = simsDF.iloc[i][key][index]

                    if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = [str(0)] * len(y)
                    elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = ["CONTROL"] * len(y)
                    else:
                        x = [str(simsDF.iloc[i][COLOR])] * len(y)

                    Y = Y + y
                    X = X + x

                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    ax.set_ylabel("VOLUME (um^3)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " VOLUME\nDISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    # ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
                    ymax = 7000
                    if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                        ymax = 700
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        plt.savefig(SAVELOC + FILEID.replace('_VOLUMES','') + '_VOLUMES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_cycles(simsDF, COLOR, FILEID, SAVELOC, TIME, HOURS):
    """Plot cell cycle distributions (in selected units) at a given time point for each cell population and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    filesplit = FILEID.split('_')

    for p in range(0, len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:
                key = 'AVG CELL CYCLES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []

                for i in range(0, len(simsDF)):
                    y = simsDF.iloc[i][key][index]

                    if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = [str(0)] * len(y)
                    elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = ["CONTROL"] * len(y)
                    else:
                        x = [str(simsDF.iloc[i][COLOR])] * len(y)

                    Y = Y + y
                    X = X + x

                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:

                    if HOURS:
                        Y = [float(c)/float(60) for c in Y]

                    figC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figC.add_subplot(1, 1, 1)
                    ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                    unit = "(hrs)" if HOURS else "(min)"
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_ylabel("CYCLE LENGTH " + unit, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " CELL CYCLE\nDISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    # ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
                    if not HOURS:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 1600
                            else:
                                ymax = 2000
                        else:
                            ymax = 2000
                    else:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 25
                            else:
                                ymax = 35
                        else:
                            ymax = 35
                    ax.set_ylim(bottom=0, top=ymax)
                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)
                    if SAVELOC == '':
                        plt.show()
                    else:
                        if HOURS:
                            units = '_HOURS'
                        else:
                            units = '_MINUTES'
                        plt.savefig(SAVELOC + FILEID.replace('_CYCLES','') + '_CYCLES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + units + '.svg', bbox_inches='tight')

    return

def plot_volumes_split(simsDF, COLOR, FILEID, SAVELOC, TIME):
    """Plot cell volume distributions at a given list of time points for each cell population and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:

                TIMES_SIM = simsDF.iloc[0]['TIME']
                index = -1

                key = 'CELL VOLUMES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []
                T = []

                for time in TIME:
                    for t in range(0, len(TIMES_SIM)):
                        if float(TIMES_SIM[t]) == float(time):
                            index = t
                            for i in range(0, len(simsDF)):
                                y = simsDF.iloc[i][key][index]
                                times = [time] * len(y)
                                if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = [str(0)] * len(y)
                                elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = ["CONTROL"] * len(y)
                                else:
                                    x = [str(simsDF.iloc[i][COLOR])] * len(y)

                                Y = Y + y
                                X = X + x
                                T = T + times

                dictVols = {'TIME': T, 'FEATURE': X, 'VOLUME': Y}
                dfVols = pd.DataFrame(dictVols, columns={'TIME': pd.Series([], dtype='float'),
                                                         'FEATURE': pd.Series([], dtype='str'),
                                                         'VOLUME': pd.Series([], dtype='float')})
                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    sns.violinplot(ax=ax, data=dfVols, x="FEATURE", y="VOLUME", order=order, hue="TIME", split=True)#, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    ax.set_ylabel("VOLUME (um^3)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " VOLUME\nDISTRIBUTIONS AT TIME " + str(TIME[0]) + " AND " + str(TIME[1]), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    ax.legend(bbox_to_anchor=(1.15, 1.0), frameon=False)
                    ymax = 7000
                    if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                        ymax = 700
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        plt.savefig(SAVELOC + FILEID.replace('_VOLUMES','') + '_VOLUMES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME[0]) + str(TIME[1]) + '.svg', bbox_inches='tight')

    return

def plot_cycles_split(simsDF, COLOR, FILEID, SAVELOC, TIME, HOURS):
    """Plot cell cycle distributions (in selected units) at a given list of time points for each cell population and color based on selected feature."""

    POP_NAMES = scripts.analyze.analyze_utilities.define_pop_names_list()
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:

                TIMES_SIM = simsDF.iloc[0]['TIME']
                index = -1

                key = 'AVG CELL CYCLES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []
                T = []

                for time in TIME:
                    for t in range(0, len(TIMES_SIM)):
                        if float(TIMES_SIM[t]) == float(time):
                            index = t
                            for i in range(0, len(simsDF)):
                                y = simsDF.iloc[i][key][index]
                                times = [time] * len(y)
                                if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = [str(0)] * len(y)
                                elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = ["CONTROL"] * len(y)
                                else:
                                    x = [str(simsDF.iloc[i][COLOR])] * len(y)

                                Y = Y + y
                                X = X + x
                                T = T + times
                if HOURS:
                    Y = [float(c) / float(60) for c in Y]

                dictVols = {'TIME': T, 'FEATURE': X, 'CYCLE': Y}
                dfVols = pd.DataFrame(dictVols, columns={'TIME': pd.Series([], dtype='float'),
                                                         'FEATURE': pd.Series([], dtype='str'),
                                                         'CYCLE': pd.Series([], dtype='float')})
                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    sns.violinplot(ax=ax, data=dfVols, x="FEATURE", y="CYCLE", order=order, hue="TIME", split=True)#, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    unit = "(hrs)" if HOURS else "(min)"
                    ax.set_ylabel("CYCLE LENGTH " + unit, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " CELL CYCLE\nDISTRIBUTIONS AT TIME " + str(TIME[0]) + " AND " + str(TIME[1]), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    ax.legend(bbox_to_anchor=(1.15, 1.0), frameon=False)

                    if not HOURS:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 1600
                            else:
                                ymax = 2000
                        else:
                            ymax = 2000
                    else:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 25
                            else:
                                ymax = 35
                        else:
                            ymax = 35
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        if HOURS:
                            units = '_HOURS'
                        else:
                            units = '_MINUTES'
                        plt.savefig(SAVELOC + FILEID.replace('_CYCLES','') + '_CYCLES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME[0]) + str(TIME[1]) + units + '.svg', bbox_inches='tight')

    return