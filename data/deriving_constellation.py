import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.table import Table

path = os.path.dirname(__file__)
datapath = os.path.dirname(path) + '/data/'


def get_cons_lines(f):
    cons_names = []
    cons_lines = []
    for line in f:
        line = line.lstrip()
        # if line.startswith(b'#'):
        #     continue
        fields = line.split()
        if not fields:
            continue
        name = fields[0]
        lines = [(int(fields[2 * i + 2]), int(fields[2 * i + 3]))
                 for i in range(int(fields[1]))]
        cons_names.append(name)
        cons_lines.append(lines)
    return cons_names, cons_lines


f = open(datapath + 'constellationship.fab', "r")

stars = Table.read(datapath + 'XHIP_7.fits')
# stars = bsc[bsc['Vmag'] < mag_lim]
stars_coord = SkyCoord(stars['RAdeg'],
                       stars['DEdeg'],
                       unit='degree',
                       frame='icrs')

cons_names, cons_lines = get_cons_lines(f)
lines = [line for lines in cons_lines for line in lines]
lines_star1 = [line[0] for line in lines]
lines_star2 = [line[1] for line in lines]
index_1 = []
for star in lines_star1:
    index_1.append(np.where(stars['HIP'] == star)[0][0])
star1_coord = stars_coord[index_1]
index_2 = []
for star in lines_star2:
    index_2.append(np.where(stars['HIP'] == star)[0][0])
star2_coord = stars_coord[index_2]
lines_xy = np.hstack((np.array([star1_coord.ra.degree,
                                star1_coord.dec.degree]).T[:, np.newaxis, :],
                      np.array([star2_coord.ra.degree,
                                star2_coord.dec.degree]).T[:, np.newaxis, :]))
print(lines_xy)
np.save(path + '/conslines.npy', lines_xy)
