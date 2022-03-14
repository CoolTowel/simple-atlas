from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, GCRS, get_moon, get_sun, get_body
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from reproject import reproject_interp
from moon import draw_moon
from astropy.visualization.wcsaxes.frame import EllipticalFrame

from astropy.visualization.wcsaxes import core

moon_radi = 1737.10 * u.km

# time
time = Time.now()
time = Time('2022-3-14 15:00:00')
# oberving location
loc_lat = 25.62
loc_lon = 101.13
loc_height = 2223
location = EarthLocation(lat=loc_lat * u.deg,
                         lon=loc_lon * u.deg,
                         height=loc_height * u.m)
obsgeoloc, obsgeovel = location.get_gcrs_posvel(time)
loc_gcrs = GCRS(obstime=time, obsgeoloc=obsgeoloc, obsgeovel=obsgeovel)

az = AltAz(obstime=time, location=location)
zenith = AltAz(az=0 * u.deg, alt=90 * u.deg, obstime=time, location=location)

# moon_eq = get_moon(time).transform_to(
#     loc_gcrs)  # moon sky coordinate in local equatorial frame
# moon_ang_radi = np.arcsin(moon_radi / moon_eq.distance.to(u.km)).to(
#     u.arcmin)  # moon angular radius

zenith_eq = zenith.transform_to(loc_gcrs)
zenith_ra = zenith_eq.ra.value
zenith_dec = zenith_eq.dec.value

FOV = 180  # field of view, diameter in degree
size_pixels = 4000
projection_scale = np.tan(FOV / 720 * np.pi) / (FOV / 720 * np.pi)
p_scale = size_pixels / projection_scale / FOV  # pixels per degree at the reference point (CD_ij(i=j) in fits)
im_w = size_pixels
im_h = size_pixels

im = np.ones((im_h, im_w))

wcs = WCS(
    header={
        'CTYPE1': 'RA---STG',
        'CTYPE2': 'DEC--STG',
        'CRVAL1': zenith_ra,  # zenith_icrs.ra.value,
        'CRPIX1': im_w / 2 + 0.5,
        'CRVAL2': zenith_dec,
        'CRPIX2': im_h / 2 + 0.5,
        'CD1_1': -1 / p_scale,
        'CD1_2': 0,
        'CD2_1': 0,
        'CD2_2': 1 / p_scale,
        'RADESYS': 'GCRS'
    })

im_r = np.full_like(im, 25 / 255)
im_g = np.full_like(im, 32 / 255)
im_b = np.full_like(im, 41 / 255)
im_a = np.full_like(im, 1)
im_rgba = np.dstack((im_r, im_g, im_b, im_a))
# im_rgba = np.dstack((im_a, im_a, im_a, im_a)) # white background

im_hdu = fits.PrimaryHDU(data=im, header=wcs.to_header())

fig = plt.figure(figsize=(5, 5), dpi=150, facecolor='#192029')
ax = fig.add_subplot(111,
                     projection=WCS(im_hdu.header, naxis=2),
                     frame_class=EllipticalFrame)

ax.imshow(im_rgba)

# moon_patch = SphericalCircle((moon_eq.ra, moon_eq.dec),
#                        moon_ang_radi * 20,
#                        edgecolor='blue',
#                        facecolor='grey',
#                        lw=0.05,
#                        ls='-',
#                        transform=ax.get_transform('world'))
horizon = SphericalCircle((zenith_ra * u.deg, zenith_dec * u.deg),
                          90 * u.deg,
                          edgecolor='grey',
                          facecolor='none',
                          lw=1,
                          ls='-',
                          transform=ax.get_transform('world'))
# ax.add_patch(horizon)
# ax.add_patch(moon_patch)

sun = get_body('sun', time, location)
moon, moon_eq = draw_moon(time, location, sun, zoom_factor=20)

moon, fo = reproject_interp(moon, im_hdu.header,order=0)
moon[np.isnan(moon)] = 0
moon_insky = np.dstack((moon,moon,moon,np.ones_like(moon)))
moon_insky[moon==0.0]=im_rgba[0,0,0:4]
moon_insky[moon==0.0,3]=0
matplot_im = ax.imshow(moon_insky)
# matplot_im.set_clip_path(ax.coords.frame.patch)
ax.text(0,90,'N',color='white', size = 10, weight ='bold',transform=ax.get_transform('world'))
ax.set_axisbelow(True)
ax.coords.frame.set_color((0.75, 0.75, 0.75))
ax.coords.frame.set_linewidth(2)
lon = ax.coords[0]
lat = ax.coords[1]
lat.set_ticks([0] * u.deg)
lat.grid(True, color='y', ls='solid',zorder=-1.0)
lon.grid(True, color='y', ls='solid',zorder=-1.0)
lon.grid(True)
lon.set_ticks_visible(False)
lon.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
lat.set_ticklabel_visible(False)
lon.set_axislabel('')
lat.set_axislabel('')

moon_azalt = moon_eq.transform_to(az)
# print((moon_azalt.az+moon_eq.position_angle(sun)-moon_eq.position_angle(zenith_eq))%360)

if __name__ == '__main__':
    plt.show()
    # ax.get_transform('world')