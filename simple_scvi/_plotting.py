import os
import numpy as np
from matplotlib import pyplot as plt

import seaborn as sns


def keys_to_colors(keys):
    uk = np.unique(keys)
    n_colors = len(uk)
    if n_colors <= 10:
        pal = sns.color_palette("tab10")[0:n_colors]
    elif n_colors <= 20:
        pal = sns.color_palette("tab20")[0:n_colors]
    else:
        raise RuntimeError("not enough colors to plot > 20 categories!")
    color_dict = dict(zip(uk, pal))
    colors = [color_dict[key] for key in keys]
    return colors


def determine_nrows_ncols(nplots: int):
    """Determine number of rows and columns a grid of subplots.
    :param nplots: total number of subplots
    :type nplots: int
    """
    if nplots < 4:
        ncols = nplots
    elif nplots < 5:
        ncols = 2
    elif nplots < 10:
        ncols = 3
    else:
        ncols = 4
    nrows = int(np.ceil(nplots / ncols))
    return nrows, ncols


def draw_plot(save_name, **kwargs):
    """Function to be used always when a plot is to be shown or saved."""
    if save_name is None:
        plt.show()
    else:
        if not os.path.exists("figures"):
            print("Creating 'figures' directory...")
            os.makedirs("figures")
        save_path = os.path.join("figures", save_name)
        plt.savefig(save_path, **kwargs)
        plt.close()


def pair_plot(
    z,
    colors=None,
    panel_size=None,
    scatter_kwargs=None,
    u=None,
    v=None,
    trajectories=None,
    traj_kwargs=None,
    quiver_kwargs=None,
    save_name=None,
    start_idx=None,
    end_idx=None,
    **kwargs
):
    """Pair plot with all dimension pairs."""
    if scatter_kwargs is None:
        scatter_kwargs = dict(alpha=0.7)
    if traj_kwargs is None:
        traj_kwargs = dict(alpha=0.7)
    if quiver_kwargs is None:
        quiver_kwargs = dict(alpha=0.7)
    d = z.shape[1]
    if d > 10:
        print("Too many pair plots! Skipping.")
        return
    nplots = int(d * (d - 1) / 2)
    if panel_size is None:
        panel_size = 6.0 if (nplots == 1) else 4.0
    nrows, ncols = determine_nrows_ncols(nplots)
    figsize = (panel_size * ncols, panel_size * nrows)
    fix, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    counter = 0
    for i in range(0, d):
        for j in range(i + 1, d):
            c = counter % ncols
            r = int(np.floor(counter / ncols))
            title = "dim " + str(i + 1) + " vs. dim " + str(j + 1)

            # Scatter plot (and possible quiver)
            if nplots > 1:
                axis = ax[r, c] if nrows > 1 else ax[c]
            else:
                axis = ax
            axis.scatter(z[:, i], z[:, j], c=colors, **scatter_kwargs)
            axis.set_title(title)
            if v is not None:
                axis.quiver(u[:, i], u[:, j], v[:, i], v[:, j], **quiver_kwargs)

            # Trajectories
            if trajectories is not None:
                n_traj = trajectories.shape[1]
                for k in range(0, n_traj):
                    zi = trajectories[:, k, i]
                    zj = trajectories[:, k, j]
                    axis.plot(zi, zj, "k", **traj_kwargs)
                    axis.plot(zi[0], zj[0], "kx", alpha=0.05)

            # Start and end index
            if start_idx is not None:
                axis.scatter(z[start_idx, i], z[start_idx, j], c="black", marker="x")
            if end_idx is not None:
                axis.scatter(z[end_idx, i], z[end_idx, j], c="red", marker="x")
            counter += 1

    # Remove extra subplots
    while counter < nrows * ncols:
        c = counter % ncols
        r = int(np.floor(counter / ncols))
        axis = ax[r, c] if nrows > 1 else ax[c]
        axis.axis("off")
        counter += 1

    draw_plot(save_name, **kwargs)


def plot_z_3d(z, colors, save_name=None, H: float = 2.5, alpha=0.8, **kwargs):
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(1, 1, 1, projection="3d")

    ax1.scatter(z[:, 0], z[:, 1], z[:, 2], c=colors, alpha=alpha)
    ax1.set_xlim(-H, H)
    ax1.set_ylim(-H, H)
    ax1.set_zlim(-H, H)

    draw_plot(save_name, **kwargs)
