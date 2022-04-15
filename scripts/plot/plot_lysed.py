import scripts.plot.plot_utilities
import matplotlib.pyplot as plt

def plot_counts_lysed_time(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot lysed cell counts over time for given population and color based on selected feature."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    if 'HEALTHY' in POP_NAME or 'DISH' in FILEID:
        xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]]
        yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
        ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0], [0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSED_' + POP_NAME.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_merge(POP_NAMES, simsDF, COLOR, FILEID, SAVELOC):
    """Plot lysed cell counts over time for both healthy and cancer cells and color based on selected feature."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    lysedLineDict = scripts.plot.plot_utilities.make_lysed_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        for POP_NAME in POP_NAMES:
            if POP_NAME != 'TISSUE':
                if 'VITRO' in FILEID:
                    plot_time = simsDF.iloc[i]['TIME']
                    lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL']
                else:
                    plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
                    lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL'][2:]
                if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict[str(0)])
                elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict["CONTROL"])
                else:
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict[str(simsDF.iloc[i][COLOR])])
    if 'HEALTHY' in POP_NAME or 'DISH' in FILEID:
        xInit = [plot_time[0], plot_time[-1]]
        yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
        ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("TISSUE CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("CELL COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_ylim(bottom=0)

    if 'VITRO' in FILEID:
        ax.set_xticks(plot_time)
    else:
        ax.set_xticks(plot_time[::10])

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for type in lysedLineDict:
        ax.plot([0], [0], color='black', label=type, linestyle=lysedLineDict[type])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSED_MERGE' + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_exact(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot lysed cell counts over time using exact death time for given population and color based on selected feature."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]*1440]
    yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
    ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (MINUTES)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])*1440])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0], [0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=90)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSEDEXACT_' + POP_NAME.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_exact_merge(POP_NAMES, simsDF, COLOR, FILEID, SAVELOC):
    """Plot lysed cell counts over time using exact death time for both healthy and cancer cells and color based on selected feature."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    lysedLineDict = scripts.plot.plot_utilities.make_lysed_line_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        for POP_NAME in POP_NAMES:
            if POP_NAME != 'TISSUE':
                    if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict[str(0)])
                    elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict["CONTROL"])
                    else:
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict[str(simsDF.iloc[i][COLOR])])
    xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]*1440]
    yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
    ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (MINUTES)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("TISSUE CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("CELL COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])*1440])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for type in lysedLineDict:
        ax.plot([0], [0], color='black', label=type, linestyle=lysedLineDict[type])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=90)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSEDEXACT_MERGE' + '.svg', bbox_inches='tight')

    return