import matplotlib.pyplot as plt 
import numpy as np


class ShowCaller:
    def show(self):
        plt.show()


def pair_plot(columns, names, units=None, groups=None, s=3, alpha=.5, groupsname=None):
    # columns are 1D iterables representing the data
    # each column should have a name and optional unit
    ncol = len(columns)
    columns = [np.array(column) for column in columns]
    group_labels = [None]
    group_indices = [np.arange(len(columns[0]))]
    if groups is not None:
        group_labels = np.unique(groups)
        group_indices = [np.where(groups == group_label) for group_label in group_labels]

    legend_handles = []
    fig, axs = plt.subplots(ncol, ncol, layout='constrained')
    for y, ycol in enumerate(columns):
        for x, xcol in enumerate(columns):
            for group_label, group_index in zip(group_labels, group_indices):
                plt.subplot(ncol, ncol, ncol*y + x+1)
                if y == x:
                    plt.hist(xcol[group_index])
                else:
                    scatter = plt.scatter(xcol[group_index], ycol[group_index], alpha=alpha, s=s)
                    if y == 0 and x == ncol-1:
                        legend_handles.append(scatter)

                if y == ncol-1 and x == 0:
                    plt.xlabel(names[x] + (f' ({units[x]})' if units else ''))
                    plt.ylabel(names[y] + (f' ({units[y]})' if units else ''))
                elif y == ncol-1:
                    plt.xlabel(names[x] + (f' ({units[x]})' if units else ''))
                    plt.yticks([])
                elif x == 0 and y > 0:
                    plt.ylabel(names[y] + (f' ({units[y]})' if units else ''))
                    plt.xticks([])
                elif y == 0 and x == ncol-1:
                    plt.ylabel(names[y] + (f' ({units[y]})' if units else ''))
                    plt.gca().yaxis.tick_right()
                    plt.gca().yaxis.set_label_position('right')
                    plt.xticks([])
                else:
                    plt.xticks([])
                    plt.yticks([])

    fig.legend(legend_handles, group_labels, loc='outside left upper', title=groupsname)
    return ShowCaller()


def matrix_bubble(M, xlabels=None, ylabels=None, cmap='RdBu', minval=None, centerval=None, maxval=None, colorbar=True, with_phase=True):
    fig, ax = plt.subplots()
    cmap = plt.get_cmap(cmap)
    nrow, ncol = M.shape
    plt.hlines(np.arange(nrow)-.5, -1, ncol, colors='grey', alpha=.5, linewidth=1)
    plt.vlines(np.arange(ncol)-.5, -1, nrow, colors='grey', alpha=.5, linewidth=1)

    minval = minval if minval is not None else M.min()
    maxval = maxval if maxval is not None else M.max()
    centerval = centerval if centerval is not None else minval

    for x in range(ncol):
        for y in range(nrow):
            v = (M[y, x] - minval)/(maxval - minval)
            c = (centerval - minval)/(maxval - minval)
            print(M[y, x], v, c)
            circle = plt.Circle((x, y), abs(v - c) * .7 + .1, color=cmap(v))
            ax.add_patch(circle)

    plt.xticks(range(ncol), xlabels, rotation=90)
    plt.yticks(range(nrow), ylabels)
    plt.xlim(-.5, ncol-.5)
    plt.ylim(-.5, nrow-.5)
    plt.gca().set_aspect('equal')
    # if colorbar:
    #   cbar = plt.colorbar(plt.gcf())
