from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np

def bv2rgb(bv): # Thanks to https://stackoverflow.com/a/43105142 
    def func(bv):
        if bv < -0.40: bv = -0.40
        if bv > 2.00: bv = 2.00

        r = 0.0
        g = 0.0
        b = 0.0

        if  -0.40 <= bv<0.00:
            t=(bv+0.40)/(0.00+0.40)
            r=0.61+(0.11*t)+(0.1*t*t)
        elif 0.00 <= bv<0.40:
            t=(bv-0.00)/(0.40-0.00)
            r=0.83+(0.17*t)
        elif 0.40 <= bv<2.10:
            t=(bv-0.40)/(2.10-0.40)
            r=1.00
        if  -0.40 <= bv<0.00:
            t=(bv+0.40)/(0.00+0.40)
            g=0.70+(0.07*t)+(0.1*t*t)
        elif 0.00 <= bv<0.40:
            t=(bv-0.00)/(0.40-0.00)
            g=0.87+(0.11*t)
        elif 0.40 <= bv<1.60:
            t=(bv-0.40)/(1.60-0.40)
            g=0.98-(0.16*t)
        elif 1.60 <= bv<2.00:
            t=(bv-1.60)/(2.00-1.60)
            g=0.82-(0.5*t*t)
        if  -0.40 <= bv<0.40:
            t=(bv+0.40)/(0.40+0.40)
            b=1.00
        elif 0.40 <= bv<1.50:
            t=(bv-0.40)/(1.50-0.40)
            b=1.00-(0.47*t)+(0.1*t*t)
        elif 1.50 <= bv<1.94:
            t=(bv-1.50)/(1.94-1.50)
            b=0.63-(0.6*t*t)

        return (r, g, b)
    r,b,g = np.vectorize(func)(np.array(bv))
    return np.vstack((r,b,g)).T

def mag2size(mag,mag_lim,max_size=10, min_size=0.05): # actually not the max size because there is no star brighter than sirius whose Vmag is -1.5
    s = max_size+(-2-np.array(mag))*(max_size-min_size)/(mag_lim+2)
    return s

def draw_star(ax, az_frame, mag_lim=5, zoom_factor=1):
    bsc = Table.read('data/Yale_Bright_Star_Catalog_V5_9095.fits')
    stars = bsc[bsc['Vmag']<mag_lim]
    stars_coord = SkyCoord(stars['_RAJ2000'], stars['_DEJ2000'], unit = 'degree', frame='icrs' ).transform_to(az_frame)
    colors = bv2rgb(stars['B-V'])
    sizes = mag2size(stars['Vmag'],mag_lim)
    ax.scatter(stars_coord.az.rad, np.pi/2-stars_coord.alt.rad, s=sizes, c=colors)