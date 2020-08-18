import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

import numpy as np
import blue2red
from Population import Population

def relabel_lam(lami_est, wi_mean):
    K = wi_mean.shape[0]
    k_ord = np.argsort(wi_mean)
    lami_new = lami_est + 0
    counts = []
    for k in range(K + 1):
        if k == 0:
            idx_k = (lami_est == 0)
        else:
            idx_k = (lami_est - 1 == k_ord[k - 1])
        lami_new[idx_k] = k
        counts.append(idx_k.sum())
    return (lami_new, counts)


def add_gridlines_Z(Z, color='grey', lw=.5):
    J, K = Z.shape
    for j in range(J):
        plt.axhline(y=j+.5, color=color, linewidth=lw)

    for k in range(K):
        plt.axvline(x=k+.5, color=color, linewidth=lw)


gridlines = add_gridlines_Z

def plot_y(yi, wi_mean, lami_est, fs_lab=10, fs_cbar=10, lw=3,
           cm=blue2red.cm(6), vlim=(-3, 3), fs_xlab=10, fs_ylab=10,
           markernames=[], interpolation=None, rotation=90, ha="center"):
    J = yi.shape[1]
    vmin, vmax = vlim

    if type(wi_mean) == int:
        K = wi_mean
        wi_mean = np.array([(lami_est == k + 1).mean() for k in range(K)])

    lami_new, counts = relabel_lam(lami_est, wi_mean)
    counts_cumsum = np.cumsum(counts)
    yi_sorted = yi[np.argsort(lami_new), :]

    im = plt.imshow(yi_sorted, aspect='auto', vmin=vmin, vmax=vmax,
                    cmap=cm, interpolation=interpolation)
    for c in counts_cumsum[:-1]:
        plt.axhline(c, color='yellow', linewidth=lw)
    plt.xticks(rotation=rotation, ha=ha)
    if len(markernames) == 0:
        plt.xticks(np.arange(J), np.arange(J) + 1, fontsize=fs_xlab)
    else:
        plt.xticks(np.arange(J), markernames, fontsize=fs_xlab)
    plt.yticks(fontsize=fs_ylab)
    plt.xlabel("markers", fontsize=fs_lab)
    plt.ylabel("cells", fontsize=fs_lab)

    ax = plt.gca()
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("top", size="7%", pad="2%")
    cax.xaxis.set_ticks_position("top")
    cbar = colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.tick_params(labelsize=fs_cbar)


def plot_Z_only(Z, fs=10, xlab=None, ylab=None, rotate_xticks=True,
                cm_greys=plt.cm.get_cmap('Greys', 5), rotation=90, ha="center"):
    plt.imshow(Z, aspect='auto', vmin=0, vmax=1, cmap=cm_greys)
    plt.xlabel(xlab, fontsize=fs)
    plt.ylabel(ylab, fontsize=fs)

    J, K = Z.shape
    plt.yticks(np.arange(J), np.arange(J) + 1, fontsize=fs)
    add_gridlines_Z(Z)
    if rotate_xticks:
        plt.xticks(rotation=rotation, fontsize=fs, ha=ha)
    else:
        plt.xticks(fontsize=fs)
    plt.xticks(np.arange(K), np.arange(K) + 1)
  

def plot_Z(Z_mean, wi_mean, lami_est, w_thresh=.01,
           cm_greys=plt.cm.get_cmap('Greys', 5), fs_lab=10,
           add_colorbar=True, fs_cbar=10, fs_w=10, fs_celltypes=10,
           xlab="markers", ylab="cell subpopulations (abundance)",
           population=0,  # type: Population
           rotation=90, ha="center",
           markernames=[], fs_markers=10, w_digits=1):

    J = Z_mean.shape[0]
    k_ord = wi_mean.argsort()
    z_cols = []

    for k in k_ord.tolist():
        if wi_mean[k] > w_thresh:
            z_cols.append(k)

    z_cols = np.array(z_cols)
    Z_hat = Z_mean[:, z_cols].T

    im = plt.imshow(Z_hat, aspect='auto', vmin=0, vmax=1, cmap=cm_greys)
    plt.xlabel(xlab, fontsize=fs_lab)
    plt.ylabel(ylab, fontsize=fs_lab)

    # W percentages
    w_perc = wi_mean[z_cols]
    w_perc = [str((wp * 100).round(w_digits)) + '%' for wp in w_perc]

    ax = plt.gca()
    # plt.xticks([])
    if population == 0:
        labels = ['{} ({})'.format(zc + 1, wp) for (zc, wp) in zip(z_cols, w_perc)]
    else:
        feature_names = [population.label(Z_mean[:, k])
                         for k in z_cols[::-1]][::-1]
        labels = ['{} ({})'.format(fname, wp)
                  for (fname, wp) in zip(feature_names, w_perc)]

    plt.yticks(np.arange(len(z_cols)), labels, fontsize=fs_celltypes)
    add_gridlines_Z(Z_hat)
    plt.xticks(rotation=rotation, fontsize=fs_markers, ha=ha)
    if len(markernames) == 0:
        plt.xticks(np.arange(J), np.arange(J) + 1)
    else:
        plt.xticks(np.arange(J), markernames)

    # add wi_mean on right side
    # K = z_cols.shape[0]
    # ax2 = ax.twinx()
    # ax2.set_yticks(range(K))
    # plt.yticks((K-1) / K * np.arange(K) + .5, w_perc[::-1], fontsize=fs_w)
    # ax2.tick_params(length=0)

    # colorbar
    if add_colorbar:
      ax_divider = make_axes_locatable(ax)
      cax = ax_divider.append_axes("top", size="7%", pad="2%")
      cax.xaxis.set_ticks_position("top")
      cbar = colorbar(im, cax=cax, orientation="horizontal")
      cbar.ax.tick_params(labelsize=fs_cbar)


def plot_yz(yi, Z_mean, wi_mean, lami_est, w_thresh=.01,
            cm_greys = plt.cm.get_cmap('Greys', 5), markernames=[],
            rotation=90, ha="center",
            cm_y=blue2red.cm(6), vlim_y=(-3, 3), fs_w=10, w_digits=1):
    J = yi.shape[1]

    vmin_y, vmax_y = vlim_y
    # cm_y.set_bad(color='black')
    # cm_y.set_under(color='blue')
    # cm_y.set_over(color='red')

    # gs = gridspec.GridSpec(1, 2, width_ratios=[2, 5]) 
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2]) 

    # Plot y
    lami_new, counts = relabel_lam(lami_est, wi_mean)
    counts_cumsum = np.cumsum(counts)
    yi_sorted = yi[np.argsort(lami_new), :]

    plt.subplot(gs[0])
    im = plt.imshow(yi_sorted, aspect='auto', vmin=vmin_y, vmax=vmax_y, cmap=cm_y)
    for c in counts_cumsum[:-1]:
        plt.axhline(c, color='yellow')
    plt.xticks(rotation=rotation, ha=ha)
    if len(markernames) == 0:
        plt.xticks(np.arange(J), np.arange(J) + 1)
    else:
        plt.xticks(np.arange(J), markernames)

    ax = plt.gca()
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("top", size="7%", pad="2%")
    cax.xaxis.set_ticks_position("top")
    colorbar(im, cax=cax, orientation="horizontal")

    # Plot Z
    k_ord = wi_mean.argsort()
    z_cols = []

    for k in k_ord.tolist():
        if wi_mean[k] > w_thresh:
            z_cols.append(k)

    z_cols = np.array(z_cols)
    Z_hat = Z_mean[:, z_cols].T

    plt.subplot(gs[1])
    im = plt.imshow(Z_hat, aspect='auto', vmin=0, vmax=1, cmap=cm_greys)
    ax = plt.gca()
    plt.xticks([])
    plt.yticks(np.arange(len(z_cols)), z_cols + 1, fontsize=fs_w)
    add_gridlines_Z(Z_hat)
    plt.colorbar(orientation='horizontal', pad=.05)

    # add wi_mean on right side
    K = z_cols.shape[0]
    ax2 = ax.twinx()
    ax2.set_yticks(range(K))
    w_perc = wi_mean[z_cols]
    w_perc = [str((wp * 100).round(w_digits)) + '%' for wp in w_perc]
    plt.yticks((K-1) / K * np.arange(K) + .5, w_perc[::-1], fontsize=fs_w)
    plt.yticks()
    ax2.tick_params(length=0)

    fig = plt.gcf()
    fig.subplots_adjust(hspace=0.2)

def colorbar_horizontal(im):
    ax = plt.gca()
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("top", size="7%", pad="2%")
    cax.xaxis.set_ticks_position("top")
    cbar = colorbar(im, cax=cax, orientation="horizontal")

def plot_y_centroids(yi, lami, wi, vlim=(-4, 4), fs_xlabel=12, fs_ylabel=12,
                     fs_xticks=12, fs_yticks=12, rotation=90, ha="center",
                     gridlines_color='black', gridlines_lw=1,
                     cm=blue2red.cm(9), population=None, Zi=None):
    J = yi.shape[1]
    K = wi.shape[0]

    selected_features = np.argwhere(wi > 0)[:, 0]
    K_sel = selected_features.shape[0]
    wi_sel = wi[selected_features]
    selected_features = selected_features[np.argsort(wi_sel)]
    wi_sel_sorted= wi_sel[np.argsort(wi_sel)]

    y_centers = np.zeros((J, K_sel))
    yticks = [0] * K_sel

    for j in range(J):
        for k in range(K_sel)[::-1]:
            k_ = selected_features[k] + 1
            y_centers[j, k] = yi[lami == k_, j].mean()
            w_ik_perc = (wi_sel_sorted[k]*100).round(1)
            if population is None or Zi is None:
                yticks[k] = '{} ({}%)'.format(k_, w_ik_perc)
            else:
                label = population.label(Zi[:, k_ - 1])
                yticks[k] = '{} ({}%)'.format(label, w_ik_perc)

    im = plt.imshow(y_centers.T, aspect='auto', cmap=cm,
                    vmin=vlim[0], vmax=vlim[1])
    plt.xticks(range(J), np.arange(J) + 1, rotation=rotation,
               fontsize=fs_xticks, ha=ha)
    plt.yticks(range(K_sel), yticks, fontsize=fs_yticks)
    plt.xlabel('markers', fontsize=fs_xlabel)
    plt.ylabel('subpopulations', fontsize=fs_ylabel)
    gridlines(y_centers.T, color=gridlines_color, lw=gridlines_lw)
    colorbar_horizontal(im)

