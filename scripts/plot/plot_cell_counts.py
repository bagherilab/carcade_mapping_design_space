import scripts.plot.plot_utilities
import matplotlib.pyplot as plt

def plot_counts(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over time for given population and color based on selected feature and indicate dose via linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over time for given population and color based on selected feature where linestyle is the same regardless of dose."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_xticks(ticks=plot_time[::2])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over time for given population and color based on selected feature and indicate live vs total populations based on linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    liveLineDict = scripts.plot.plot_utilities.make_live_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 100000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_treat(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over treatment time for given population and color based on selected feature and indicate dose via linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_treat_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over treatment time for given population and color based on selected feature where linestyle is the same regardless of dose."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_treat_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts over treatment time for given population and color based on selected feature and indicate live vs total populations based on linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    liveLineDict = scripts.plot.plot_utilities.make_live_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE'][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    liveLineDict = scripts.plot.plot_utilities.make_live_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            killed_live = [simsDF.iloc[i][POP_NAME + ' LIVE'][a] / simsDF.iloc[i][POP_NAME + ' LIVE'][0] for a in
                      range(0, len(simsDF.iloc[i][POP_NAME + ' LIVE']))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            killed_live = [simsDF.iloc[i][POP_NAME + ' LIVE'][a] / simsDF.iloc[i][POP_NAME + ' LIVE'][44] for a in
                      range(44, len(simsDF.iloc[i][POP_NAME + ' LIVE']))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_norm(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts normalized by count/dose at start of treatment over time for given population and color based on selected feature and indicate dose via linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:
            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1 - frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
            else:
                plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_norm_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts normalized by count/dose at start of treatment over time for given population and color based on selected feature where linestyle is the same regardless of dose."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]
    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:

            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1.0-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
            else:
                plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_norm_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot cell counts normalized by count/dose at start of treatment over time for given population and color based on selected feature and indicate live vs total populations based on linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    liveLineDict = scripts.plot.plot_utilities.make_live_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:
            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']

                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
                elif POP_NAME in ['T-CELL']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/simsDF.iloc[i][POP_NAME + ' LIVE'][0] for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
            else:
                plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD8']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]
                elif POP_NAME in ['T-CELL']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/simsDF.iloc[i][POP_NAME + ' LIVE'][2] for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict[str(0)])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict["CONTROL"])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_CH_scatter(simsDF, COLOR, MARKER, FILEID, SAVELOC, TIME):
    """Plot scatter plot of final living cancer vs healthy cell counts per simulation."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    MARKER_DICT = scripts.plot.plot_utilities.make_features_marker_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]
    markerDict = MARKER_DICT[MARKER]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCH = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCH.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        cancer = simsDF.iloc[i]['CANCER LIVE'][index]
        healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]
        if int(simsDF.iloc[i]['DOSE']) == 0:
            color = 'black'
        else:
            color = colorDict[str(simsDF.iloc[i][COLOR])]
        if MARKER != 'o':
            marker = markerDict[simsDF.iloc[i][MARKER]]
        else:
            marker = markerDict[MARKER]

        ax.scatter(cancer, healthy, color=color, marker=marker, s=150) # facecolor='none', edgecolor=color, marker=marker, s=200, linewidth=2)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    if MARKER != 'o':
        for marker in markerDict:
            ax.plot([0],[0], color='black', label=marker, linestyle=markerDict[marker])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)
    ax.set_title("LIVING CANCER VS HEALTHY CELLS", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("LIVING CANCER CELL COUNT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("LIVING HEALTHY CELL COUNT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_SCATTER_CH_' + COLOR.replace(' ','') + '_' + MARKER.replace(' ','') + '_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_CH_scatter_normalized(simsDF, COLOR, MARKER, FILEID, SAVELOC, TIME):
    """Plot scatter plot of final living cancer vs healthy cell counts normalized to start of treatment per simulation."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    MARKER_DICT = scripts.plot.plot_utilities.make_features_marker_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    colorDict = COLOR_DICT[COLOR]
    markerDict = MARKER_DICT[MARKER]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCH = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCH.add_subplot(1, 1, 1)

    if 'VITRO' in FILEID:
        xmax = 6
        if COLOR == 'CAR AFFINITY' and '_CH_' in FILEID:
            ymax = 1.4
        elif COLOR == 'ANTIGENS CANCER' and '_CH_' in FILEID:
            ymax = 1.4
        else:
            ymax = 1.4
    else:
        xmax = 2.5
        ymax = 1


    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            cancer = simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][0]
            healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]/simsDF.iloc[i]['HEALTHY LIVE'][0]
        else:
            cancer = simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][44]
            healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]/simsDF.iloc[i]['HEALTHY LIVE'][44]
        if int(simsDF.iloc[i]['DOSE']) == 0:
            color = 'black'
        else:
            color = colorDict[str(simsDF.iloc[i][COLOR])]
        if MARKER != 'o':
            marker = markerDict[simsDF.iloc[i][MARKER]]
        else:
            marker = markerDict[MARKER]

        ax.scatter(cancer, healthy,  color=color, marker=marker, s=150) # facecolor='none', edgecolor=color, marker=marker, s=200, linewidth=2)

    ax.plot([1, 1], [0, ymax], color='lightgray', zorder=-1, linewidth=5)
    ax.plot([0, xmax], [1, 1], color='lightgray', zorder=-1, linewidth=5)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    if MARKER != 'o':
        for marker in markerDict:
            ax.plot([0],[0], color='black', label=marker, linestyle=markerDict[marker])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)
    ax.set_title("LIVING CANCER VS HEALTHY\nCELLS NORMALIZED", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("LIVING CANCER CELL\nCOUNT NORMALIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("LIVING HEALTHY CELL\nCOUNT NORMALIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0,xmax])
    ax.set_ylim([0,ymax])

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_SCATTERNORM_CH_' + COLOR.replace(' ','') + '_' + MARKER.replace(' ','') + '_' + str(TIME) + '.svg', bbox_inches='tight')

    return