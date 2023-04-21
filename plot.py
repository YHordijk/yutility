import matplotlib.pyplot as plt 
import matplotlib
from yutility import ensure_list
import numpy as np
import atexit


class ShowCaller:
    def show(self, block=False):
        atexit.register(plt.show, block=True)
        plt.show(block=block)
        return self

    def savefig(self, *args, **kwargs):
        plt.savefig(*args, **kwargs)
        return self

    def close(self):
        plt.close()


def density(points, lim, resolution=1000, s2=.002):
    s2 = s2 * (lim[1] - lim[0])**2
    x = np.linspace(*lim, resolution)
    dens = np.zeros(resolution)
    for point in points:
        dens += np.exp(-(x - point)**2/s2)
    return x, dens/dens.max()


def scatter(x, y, xlabel=None, ylabel=None, plot_marginals=True, s=3, alpha=.5, groups=None, groupsname=None, legendloc='outside right upper'):
    # Create Fig and gridspec
    if plot_marginals:
        grid = plt.GridSpec(4, 4, hspace=0, wspace=0)
    else:
        grid = plt.GridSpec(1, 1)

    group_labels = [None]
    group_indices = [np.arange(len(x))]
    if groups is not None:
        group_labels = np.unique(groups)
        group_indices = [np.where(groups == group_label) for group_label in group_labels]

    # Define the axes
    if plot_marginals:
        ax_main = plt.gcf().add_subplot(grid[1:, :-1])
        ax_right = plt.gcf().add_subplot(grid[1:, -1], xticks=[], yticks=[], xticklabels=[], yticklabels=[])
        ax_top = plt.gcf().add_subplot(grid[0, :-1], xticks=[], yticks=[], xticklabels=[], yticklabels=[])
    else:
        ax_main = plt.gcf().add_subplot(grid[0, 0])

    for group_label, group_index in zip(group_labels, group_indices):
        ax_main.scatter(x[group_index], y[group_index], s=s, alpha=alpha, label=group_label)

    if plot_marginals:
        for group_label, group_index in zip(group_labels, group_indices):
            ylim = ax_main.get_ylim()
            ymarg = (max(ylim) - min(ylim)) * ax_main.margins()[1]
            ylim = min(ylim) - ymarg, max(ylim) + ymarg

            xlim = ax_main.get_xlim()
            xmarg = (max(xlim) - min(xlim)) * ax_main.margins()[0]
            xlim = min(xlim) - xmarg, max(xlim) + xmarg

            dens = density(y[group_index], ylim)
            ax_right.plot(*dens[::-1], linewidth=1)
            ax_right.fill_betweenx(dens[0], dens[1], 0, alpha=.3)
            ax_right.margins(0)
            ax_right.axis('off')
            lim = ax_right.get_xlim()
            lim = [lim[0], lim[1] * 1.05]
            ax_right.set_xlim(lim)

            dens = density(x[group_index], xlim)
            ax_top.plot(*dens, linewidth=1)
            ax_top.fill_between(dens[0], 0, dens[1], alpha=.3)
            ax_top.margins(0)
            ax_top.axis('off')
            lim = ax_top.get_ylim()
            lim = [lim[0], lim[1] * 1.05]
            ax_top.set_ylim(lim)

    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)

    plt.gcf().legend(loc=legendloc, title=groupsname)

    return ShowCaller()


def heatmap(M, extent=None, xlabel=None, ylabel=None, plot_marginals=True, cmap='Blues', groupcmap='tab10'):
    Ms = ensure_list(M)
    # Create Fig and gridspec
    if plot_marginals:
        grid = plt.GridSpec(4, 4, hspace=0, wspace=0)
    else:
        grid = plt.GridSpec(1, 1)

    # Define the axes
    if plot_marginals:
        ax_main = plt.gcf().add_subplot(grid[1:, :-1])
        ax_right = plt.gcf().add_subplot(grid[1:, -1], xticks=[], yticks=[], xticklabels=[], yticklabels=[])
        ax_top = plt.gcf().add_subplot(grid[0, :-1], xticks=[], yticks=[], xticklabels=[], yticklabels=[])
    else:
        ax_main = plt.gcf().add_subplot(grid[0, 0])

    for i, M in enumerate(Ms):
        if plot_marginals:
            ax_main.margins(0)

            margin_right = M.sum(axis=1)
            margin_right = (margin_right - margin_right.min())/(margin_right.max() - margin_right.min())
            ax_right.plot(margin_right, range(len(margin_right)), linewidth=1)
            ax_right.fill_betweenx(range(len(margin_right)), margin_right, 0, alpha=.3)
            ax_right.margins(0)
            ax_right.axis('off')
            lim = ax_right.get_xlim()
            ax_right.set_xlim([lim[0], 1.05])

            margin_top = M.sum(axis=0)
            margin_top = (margin_top - margin_top.min())/(margin_top.max() - margin_top.min())
            ax_top.plot(range(len(margin_top)), margin_top, linewidth=1)
            ax_top.fill_between(range(len(margin_top)), margin_top, 0, alpha=.3)
            ax_top.margins(0)
            ax_top.axis('off')
            lim = ax_top.get_xlim()
            ax_top.set_ylim([lim[0], 1.05])

        imcmap = matplotlib.colors.LinearSegmentedColormap.from_list('', [(1, 1, 1, 1), plt.get_cmap(groupcmap)(i)])
        ax_main.imshow(M, origin='lower', aspect='auto', cmap=imcmap, extent=extent, alpha=M)

    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)

    return ShowCaller()


def pair_plot(columns, names, units=None, groups=None, s=3, alpha=.5, groupsname=None, figsize=None):
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
    fig, axs = plt.subplots(ncol, ncol, figsize=figsize)
    plt.subplots_adjust(hspace=0.01, wspace=0.01)
    for y, ycol in enumerate(columns):
        for x, xcol in enumerate(columns):
            for group_label, group_index in zip(group_labels, group_indices):
                if x == y:
                    continue
                ax = plt.subplot(ncol, ncol, ncol*y + x + 1)
                ax.margins(0.1)
                scatter = ax.scatter(xcol[group_index], ycol[group_index], alpha=alpha, s=s)
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

    for x, xcol in enumerate(columns):
        for group_label, group_index in zip(group_labels, group_indices):
            ax = plt.subplot(ncol, ncol, ncol*x + x + 1)
            xlim = xcol.max(), xcol.min()
            dens = density(xcol[group_index], xlim)
            ax.plot(*dens, linewidth=1)
            ax.fill_between(dens[0], 0, dens[1], alpha=.3)
            if x == ncol - 1:
                plt.yticks([])
                plt.xlabel(names[x] + (f' ({units[x]})' if units else ''))
            else:
                plt.xticks([])
                plt.yticks([])
            lim = ax.get_xlim()
            ax.set_ylim([0, 1.05])
            ax.margins(0)

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
