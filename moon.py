from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, GCRS, get_moon, get_sun, get_body
from astropy.time import Time
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
# from sympy import symbols
# from sympy.solvers import solve
from scipy import ndimage
from reproject import reproject_interp

def draw_moon(time, location, sun, zoom_factor = 1, light_color = 0.75, shadow_color = 0.25, size=1001):
    moon = get_moon(time,location)
    moon_radi = (1737.10*u.km).to(u.au)
    moon_ang_radi = np.arcsin(moon_radi/moon.distance.to(u.au)).to(u.deg) # moon's angular radius in degree
    
    FOV = moon_ang_radi.value*2*zoom_factor # field of view (also the enlarged moon size, if zoom_factor is not 1), diameter in degree (but we will add some blank pixels around it, so it is not the final FOV)
    moon_size_pixels = size # moon diameter in fits image 
    projection_scale = np.tan(FOV/720*np.pi)/(FOV/720*np.pi) # projection distortion compensation
    p_scale = moon_size_pixels/projection_scale/FOV # pixels per degree at the reference point (CD_ij(i=j) in fits)
    moon_w = moon_size_pixels
    moon_h = moon_size_pixels

    moon_radius_pixels = moon_size_pixels/2 
    terminator_radius_pixels = moon_radius_pixels*1.005 # making the radius of terminator slighly larger than moon angular size to make simulate the true moon phase

    cen = int(moon_size_pixels/2) # image center
    Y, X = np.mgrid[0:moon_h, 0:moon_w]-cen # make a larger image box
    moon_phase = 180*u.deg-np.arccos((moon.distance**2+moon.separation_3d(sun)**2-sun.distance**2)/(2*moon.distance*moon.separation_3d(sun))).to(u.deg) # 0-180deg, 0 is the mew moon, and 180 is the full moon
    moon_light = np.logical_and(X**2+Y**2<moon_radius_pixels**2,Y>(np.sqrt(terminator_radius_pixels**2-X**2)*np.cos(moon_phase))) # moon light side mask 
    moon_shadow = np.logical_and(X**2+Y**2<moon_radius_pixels**2,Y<(np.sqrt(terminator_radius_pixels**2-X**2)*np.cos(moon_phase))) # moon shadow side mask 
    moon_image = np.zeros((moon_h, moon_w))

    moon_image[moon_light] = light_color #[light_color,light_color,light_color,1]
    moon_image[moon_shadow] = shadow_color #[shadow_color,shadow_color,shadow_color,1]
    moon_image = ndimage.rotate(moon_image,360-moon.position_angle(sun).degree, reshape=False,order=0)
    pad=5
    box_w = moon_w+2*pad
    box_h = moon_h+2*pad
    box = np.zeros((box_h,box_w))
    box[pad:pad+moon_h,pad:pad+moon_w] = moon_image

    box_wcs = WCS(header={
        'CTYPE1': 'RA---STG',
        'CTYPE2': 'DEC--STG',
        'CRVAL1': moon.ra.value,
        'CRPIX1': box_w/2+0.5,
        'CRVAL2': moon.dec.value,
        'CRPIX2': box_h/2+0.5,
        'CD1_1': -1/p_scale,
        'CD1_2': 0,
        'CD2_1': 0,
        'CD2_2': 1/p_scale,
        'RADESYS': 'GCRS'
    })
    hdu = fits.PrimaryHDU(data = box, header = box_wcs.to_header())
    return hdu, moon

if __name__ == '__main__':
    time = Time('2022-3-11 14:51:21')
    # oberving location 
    loc_lat = 25.62
    loc_lon = 101.13
    loc_height = 2223
    location = EarthLocation(lat=loc_lat*u.deg, lon=loc_lon*u.deg, height=loc_height*u.m)
    obsgeoloc, obsgeovel = location.get_gcrs_posvel(time)
    loc_gcrs = GCRS(obstime=time, obsgeoloc=obsgeoloc, obsgeovel=obsgeovel)
    sun = get_body('sun', time, location)
    moon, moon_eq = draw_moon(time, location, sun)
    fig = plt.figure(figsize = (5,5), dpi=100)
    ax=fig.add_subplot(111, projection=WCS(moon.header))
    ax.imshow(moon.data,cmap='gray')
    star = SkyCoord.from_name('hip29390').transform_to(loc_gcrs)
    star = SphericalCircle((star.ra, star.dec), 0.5*u.arcmin, zorder=50, edgecolor='blue', facecolor='blue', lw=1, ls='-', transform=ax.get_transform('world'))
    ax.add_patch(star)
    moon_phase = np.arccos((moon_eq.distance**2+moon_eq.separation_3d(sun)**2-sun.distance**2)/(2*moon_eq.distance*moon_eq.separation_3d(sun))).to(u.deg)
    plt.show()
    print(moon_phase)