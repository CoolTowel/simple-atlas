from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import os

path = os.path.dirname(__file__)
datapath = os.path.dirname(path) + '/data/'


def bv2rgb(stars, key='B-V'):  # Thanks to https://stackoverflow.com/a/43105142

    def func(bv):
        if bv < -0.40: bv = -0.40
        if bv > 1.98: bv = 1.98

        r = 0.0
        g = 0.0
        b = 0.0

        if -0.40 <= bv < 0.00:
            t = (bv + 0.40) / (0.00 + 0.40)
            r = 0.61 + (0.11 * t) + (0.1 * t * t)
        elif 0.00 <= bv < 0.40:
            t = (bv - 0.00) / (0.40 - 0.00)
            r = 0.83 + (0.17 * t)
        elif 0.40 <= bv < 2.10:
            t = (bv - 0.40) / (2.10 - 0.40)
            r = 1.00
        if -0.40 <= bv < 0.00:
            t = (bv + 0.40) / (0.00 + 0.40)
            g = 0.70 + (0.07 * t) + (0.1 * t * t)
        elif 0.00 <= bv < 0.40:
            t = (bv - 0.00) / (0.40 - 0.00)
            g = 0.87 + (0.11 * t)
        elif 0.40 <= bv < 1.60:
            t = (bv - 0.40) / (1.60 - 0.40)
            g = 0.98 - (0.16 * t)
        elif 1.60 <= bv < 2.00:
            t = (bv - 1.60) / (2.00 - 1.60)
            g = 0.82 - (0.5 * t * t)
        if -0.40 <= bv < 0.40:
            t = (bv + 0.40) / (0.40 + 0.40)
            b = 1.00
        elif 0.40 <= bv < 1.50:
            t = (bv - 0.40) / (1.50 - 0.40)
            b = 1.00 - (0.47 * t) + (0.1 * t * t)
        elif 1.50 <= bv < 1.94:
            t = (bv - 1.50) / (1.94 - 1.50)
            b = 0.63 - (0.6 * t * t)

        return (r, g, b)

    bv = np.array(stars[key])
    r, b, g = np.vectorize(func)(np.array(bv))
    return np.vstack((r, b, g)).T


def mag2size(stars, key='Vmag', k_size=0.4, max_size=40, min_size=0.3):
    ''' A non-linear convert vmag to size, when vmag = -2, size = max_szie. In fact the max size would never be reached because there is no star brighter than sirius whose Vmag is -1.5.'''
    mag = np.array(stars[key])
    mag = mag - mag.min()  # normalize
    mag_lim = mag.max()
    a = (min_size - max_size) / (np.exp(-k_size * mag_lim) - 1)
    s = a * (np.exp(-k_size * mag) - 1) + max_size
    return s


def draw_star(ax,
              az_frame,
              mag_lim=5,
              k_size=0.4,
              zoom_factor=1,
              constline=True):
    stars = Table.read(datapath + 'XHIP_7.fits')
    # stars = bsc[bsc['Vmag'] < mag_lim]
    stars_coord = SkyCoord(stars['RAdeg'],
                           stars['DEdeg'],
                           unit='degree',
                           frame='icrs').transform_to(az_frame)
    mask = stars['Vmag'] < mag_lim
    colors = bv2rgb(stars[mask])
    sizes = mag2size(stars[mask], k_size=k_size) * zoom_factor
    ax.scatter(stars_coord.az.rad[mask],
               np.pi / 2 - stars_coord.alt.rad[mask],
               s=sizes,
               c=colors,
               linewidths=0.0,
               edgecolors=None,
               zorder=500)
