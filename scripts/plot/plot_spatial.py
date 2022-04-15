import scripts.plot.plot_utilities
import matplotlib.pyplot as plt

def plot_counts_radius(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):
    """Plot cell counts (normalized to locations at radius) at each radius for given cell population and color based on selected feature and indicate dose via linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) +  '.svg', bbox_inches='tight')

    return

def plot_counts_radius_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):
    """Plot cell counts (normalized to locations at radius) at each radius for given cell population and color based on selected feature where linestyle is the same regardless of dose."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) +  '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_radius_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):
    """Plot cell counts (normalized to locations at radius) at each radius for given cell population and color based on selected feature and indicate live vs total populations based on linestyle."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    liveLineDict = scripts.plot.plot_utilities.make_live_line_dict()

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_MERGE_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return