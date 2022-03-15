from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, GCRS, get_moon, get_body
from astropy.time import Time
from astropy.wcs import WCS
from astropy.table import Table
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
import matplotlib.font_manager as fm
from moon import draw_moon
from moon_logo import draw_moon_logo
from stereographic import StereographicZenithScale
import time as pytime
start_time = pytime.time()


font_path = './fonts/Libertinus-7.040/static/OTF/'
font_files = fm.findSystemFonts(fontpaths=font_path)
for font_file in font_files:
    fm.fontManager.addfont(font_file)
    
mscale.register_scale(StereographicZenithScale)
red = '#C80815'
light_blue = '#4673bc'
dark_yellow = '#d4b300'

def skychart(time=None, location=None, show=False):
    '''
    Return skychart at a time in a location. Time should be an astropy.time.Time instance and location should be an astropy.coordinates.EarthLocation instance.'''
    # time
    if time is None:
        time = Time.now()
        # time = Time('2022-2-28 00:00:00')
    # oberving location
    if location is None:
        loc_lat = 25.62
        loc_lon = 101.13
        loc_height = 2223
        location = EarthLocation(lat=loc_lat * u.deg,
                            lon=loc_lon * u.deg,
                            height=loc_height * u.m)
    
    obsgeoloc, obsgeovel = location.get_gcrs_posvel(time)
    loc_gcrs = GCRS(obstime=time, obsgeoloc=obsgeoloc, obsgeovel=obsgeovel)

    az_frame = AltAz(obstime=time, location=location)
    # zenith = AltAz(az=0 * u.deg, alt=90 * u.deg, obstime=time, location=location)

    # the celestial body
    sun = get_body('sun', time, location)
    moon = get_moon(time, location)

    fig = plt.figure(figsize=(8,8), dpi=100)
    ax = fig.add_subplot(111, projection='polar')
    ax.set_rscale('stereographiczenith')
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_yticks([])
    ax.set_ylim(0,np.pi/2)
    ax.set_theta_zero_location("N")
    ax.set_xticks([0, np.pi/2,np.pi,np.pi/2*3], ('N', 'E', 'S', 'W'))
    ax.get_xticklabels()[0].set_color(red)

    # zenith
    # the marker '+' looks ugly in polar projection
    # ax.plot(0,0, marker='+')
    ax.plot([0,np.pi],[0.04,0.04], c=dark_yellow, lw=0.5,zorder=0)
    ax.plot([np.pi/2,np.pi/2*3],[0.04,0.04], c=dark_yellow, lw=0.5,zorder=0)

    # equator
    eq_lon=np.arange(0,360,1)
    eq_lat = np.zeros_like(eq_lon)
    eq = SkyCoord(eq_lon,eq_lat,unit=u.deg,frame=loc_gcrs).transform_to(az_frame)
    eq = np.array([eq.az.value, 90 - eq.alt.value])/180 * np.pi
    ax.plot(eq[0],eq[1],zorder=-1,c=light_blue,lw=plt.rcParams['axes.linewidth'])

    # moon
    l_p, s_p = draw_moon_logo(az_frame, moon, sun, loc_gcrs,zoom_factor=15, resolution=40, light_color=0.9)
    ax.add_patch(l_p)
    ax.add_patch(s_p)

    # Bright star
    bsc = Table.read('data/Yale_Bright_Star_Catalog_V5_9095.fits')
    stars_3 = bsc[bsc['Vmag']<3]
    stars = SkyCoord(stars_3['_RAJ2000'], stars_3['_DEJ2000'], unit = 'degree', frame='icrs' ).transform_to(az_frame)
    ax.scatter(stars.az.rad, np.pi/2-stars.alt.rad,s=0.1, c='white')

    plt.tight_layout()
    if show is True:
        plt.show()
    
    return fig

if __name__ == '__main__':
    fig = skychart(show= True)
    save_path = './test/test.png'
    # fig.savefig(save_path,dpi=300)
    print("--- %s seconds ---" % (pytime.time() - start_time))