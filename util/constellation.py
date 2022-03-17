import numpy as np
import os
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection

path = os.path.dirname(__file__)
datapath = os.path.dirname(path) + '/data/'


def draw_constellation(
    ax,
    az_frame,
    lineswitdhs=0.6,
    color=[0.4] * 3,
):
    lines = np.load(datapath + 'conslines.npy')
    star1_coord = SkyCoord(lines[:, 0, 0], lines[:, 0, 1],
                           unit='degree').transform_to(az_frame)
    star2_coord = SkyCoord(lines[:, 1, 0], lines[:, 1, 1],
                           unit='degree').transform_to(az_frame)
    lines_xy = np.hstack(
        (np.array([star1_coord.az.rad,
                   np.pi / 2 - star1_coord.alt.rad]).T[:, np.newaxis, :],
         np.array([star2_coord.az.rad,
                   np.pi / 2 - star2_coord.alt.rad]).T[:, np.newaxis, :]))
    ax.add_collection(
        LineCollection(lines_xy,
                       linewidths=lineswitdhs,
                       colors=color,
                       zorder=499))
