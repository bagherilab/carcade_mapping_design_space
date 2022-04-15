import scripts.plot.plot_utilities
import matplotlib.pyplot as plt

def plot_rank_parity(rankDF, COLOR, XRANK, SAVELOC):
    """Make parity plot of rank in dish vs tissue for identical simulation setups."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)
    ax.plot([0, len(rankDF)+1], [0, len(rankDF)+1], color='lightgray', zorder=0)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='RANK DISH', ascending=True)
        x = rankDFsorted['RANK DISH']
        y = rankDFsorted['RANK TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='RANK TISSUE', ascending=True)
        x = rankDFsorted['RANK TISSUE']
        y = rankDFsorted['RANK DISH']
        for i in range(0, len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("RANK IN " + XRANK + " vs RANK IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("RANK IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("RANK IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0,len(rankDF)+1])
    ax.set_ylim([0,len(rankDF)+1])

    plt.xticks([i for i in range(1,15)], fontsize=TICKSIZE-10, rotation=45)
    plt.yticks([i for i in range(1,15)], fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'RANK_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_rank_ladder(rankDF, COLOR, XRANK, SAVELOC):
    """Make ladder plot of rank in dish and tissue for identical simulation setups."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='RANK DISH', ascending=False)
        x = rankDFsorted['RANK DISH']
        y = rankDFsorted['RANK TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1, markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='RANK TISSUE', ascending=False)
        x = rankDFsorted['RANK TISSUE']
        y = rankDFsorted['RANK DISH']
        for i in range(0, len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1, markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("RANK IN " + XRANK + " vs RANK IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("CONTEXT", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("RANK", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0.7,2.3])
    ax.set_ylim([0,len(rankDF)+1])

    plt.xticks([1, 2], labels=[XRANK, YRANK], fontsize=TICKSIZE-10)
    #plt.xlabel([XRANK, YRANK], fontsize=TICKSIZE-10, rotation=45)
    plt.yticks([i for i in range(1,15)], fontsize=TICKSIZE-10)
    plt.gca().invert_yaxis()

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'RANK_LADDER_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_parity(rankDF, COLOR, XRANK, SAVELOC):
    """Make parity plot of score in dish vs tissue for identical simulation setups."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        x = rankDFsorted['SCORE DISH']
        y = rankDFsorted['SCORE TISSUE']
        ax.plot([0, 2], [0, 2], color='lightgray', zorder=0)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 2])
        for i in range(0,len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        x = rankDFsorted['SCORE TISSUE']
        y = rankDFsorted['SCORE DISH']
        ax.plot([0, 2], [0, 2], color='lightgray', zorder=0)
        ax.set_xlim([0, 2])
        ax.set_ylim([0, 2])
        for i in range(0, len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("SCORE IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE-10)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORE_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_norm_parity(rankDF, COLOR, XRANK, SAVELOC):
    """Make parity plot of min-max normalized score values in dish vs tissue for identical simulation setups."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)
    ax.plot([0, 1], [0, 1], color='lightgray', zorder=0)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        a = rankDFsorted['SCORE DISH']
        b = rankDFsorted['SCORE TISSUE']
        x = [(i - min(a)) / (max(a) - min(a)) for i in a]
        y = [(i - min(b)) / (max(b) - min(b)) for i in b]
        for i in range(0,len(x)):
            ax.scatter(x[i], y[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        a = rankDFsorted['SCORE TISSUE']
        b = rankDFsorted['SCORE DISH']

        x = [(i - min(a)) / (max(a) - min(a)) for i in a]
        y = [(i - min(b)) / (max(b) - min(b)) for i in b]
        for i in range(0, len(x)):
            ax.scatter(x[i], y[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("SCORE IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE-10, rotation=45)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORENORM_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_ladder(rankDF, COLOR, XRANK, SAVELOC):
    """Make ladder plot of score in dish vs tissue for identical simulation setups."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        x = rankDFsorted['SCORE DISH']
        y = rankDFsorted['SCORE TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1,
                    markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        x = rankDFsorted['SCORE TISSUE']
        y = rankDFsorted['SCORE DISH']
        for i in range(0, len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1,
                    markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("CONTEXT", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0.7, 2.3])
    ax.set_ylim([0, 2])
    plt.yticks([0, 0.5, 1, 1.5, 2])
    plt.xticks([1, 2], labels=[XRANK, YRANK], fontsize=TICKSIZE - 10)
    plt.xticks(fontsize=TICKSIZE-10)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORE_LADDER_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return