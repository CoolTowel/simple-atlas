import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, position_angle, get_moon 
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import os

path = os.path.dirname(__file__)
datapath = os.path.dirname(path)+'/data/'

moon_radi = (1737.10 * u.km).to(u.au)


def moon_edge(resolution=50):
    '''
    The outline of moon. return an array of its pixel coordinates [X,Y] with 2 column
    '''
    t = np.linspace(0, np.pi * 2.0, resolution)
    t = t.reshape((len(t), 1))
    x = np.cos(t)
    y = np.sin(t)
    return np.hstack((x, y))


def moon_terminator(resolution=50, b=0.5, k=1.005):
    '''
    The ellipse of moon's terminator. A factor k was applying to the radius to make the final looks more real. return an array of its pixel coordinates [X,Y] with 2 column. The points are not uniform in the arc. 
    '''
    t = np.linspace(0, np.pi * 2.0, resolution)
    t = t.reshape((len(t), 1))
    x = np.cos(t) * k
    y = np.sin(t) / b * k
    return np.hstack((x, y))


def draw_moon_logo(ax,
                   az_frame,
                   moon,
                   sun,
                   loc_gcrs,
                   lms,
                   zoom_factor=20,
                   resolution=50,
                   light_color=np.array([1.0, 253 / 255, 230 / 255]),
                   shadow_color=0.25,
                   k=1.005):
    '''
    Get the matplotlib patch of moon's light and shadow in the latitude-zenith coordinate. Return in tuple
    '''
    moon_lat = moon.transform_to(az_frame).alt.degree
    if moon_lat<-10:
        return None

    moon_ang_radi = np.arcsin(moon_radi / moon.distance.to(u.au)).to(u.deg)

    FOV = moon_ang_radi.value * 2 * zoom_factor  # field of view for our imaginary fits image(also the enlarged moon size, if zoom_factor is not 1), diameter in degree
    moon_size_pixels = 2  # conviniant for make circle of moon_edge, so that the circle radius will be 1. do not change
    projection_scale = np.tan(FOV / 720 * np.pi) / (
        FOV / 720 * np.pi)  # projection distortion compensation
    p_scale = moon_size_pixels / projection_scale / FOV  # pixels per degree at the reference point (CD_ij(i=j) in fits)

    moon_wcs = WCS(
        header={
            'CTYPE1': 'RA---STG',
            'CTYPE2': 'DEC--STG',
            'CRVAL1': moon.ra.value,
            'CRPIX1': 1,
            'CRVAL2': moon.dec.value,
            'CRPIX2': 1,
            'CD1_1': -1 / p_scale,
            'CD1_2': 0,
            'CD2_1': 0,
            'CD2_2': 1 / p_scale,
            'RADESYS':
            'ICRS',  # actually GCRS but astropy cannot convert pixels' coordinates to world coordinates. So just use the fake one, pretending we are in ICRS.
        })

    moon_phase = 180 * u.deg - np.arccos(
        (moon.distance**2 + moon.separation_3d(sun)**2 - sun.distance**2) /
        (2 * moon.distance * moon.separation_3d(sun))).to(u.deg)
    # 0-180deg, 0 is the mew moon, and 180 is the full moon
    light_side_angle = moon.position_angle(sun).rad
    b = np.absolute(1 / (np.cos(moon_phase)))
    edge = moon_edge(resolution=resolution)
    termi = moon_terminator(resolution=resolution, k=k, b=b)
    # k, terminator radius scale. Making the radius of terminator slighly larger than moon angular size to make simulate the true moon phase
    cross_y = np.sqrt((k**2 - 1) / (b**2 - 1))

    if np.cos(moon_phase) > 0:
        if cross_y >= 1:
            light = None
            shadow = edge
        else:
            cross_x = np.sqrt((b**2 - k**2) / (b**2 - 1))
            moon_edge_light = edge[edge[:, 1] > cross_y, :]
            moon_edge_light = np.vstack(([cross_x, cross_y], moon_edge_light))
            moon_terminator_line = termi[termi[:, 1] > cross_y, :]
            moon_terminator_line = np.vstack(
                ([-cross_x, cross_y], moon_terminator_line[::-1]))
            light = np.vstack((moon_edge_light, moon_terminator_line))
            indices = np.where(edge[:, 1] >= cross_y)[0]
            lower = 1 - (len(edge[:, 0]) - indices[-1])
            upper = indices[0]
            moon_edge_shadow = np.vstack((edge[lower:, :], edge[:upper, :]))
            moon_edge_shadow = np.vstack(([cross_x,
                                           cross_y], moon_edge_shadow[::-1]))
            shadow = np.vstack((moon_terminator_line, moon_edge_shadow))
    elif np.cos(moon_phase) < 0:
        if cross_y >= 1:
            light = edge
            shadow = None
        else:
            cross_x = np.sqrt((b**2 - k**2) / (b**2 - 1))
            cross_y = -cross_y
            indices = np.where(edge[:, 1] <= cross_y)[0]
            lower = 1 - (len(edge[:, 0]) - indices[-1])
            upper = indices[0]
            moon_edge_light = np.vstack((edge[lower:, :], edge[:upper, :]))
            moon_edge_light = np.vstack(([cross_x, cross_y], moon_edge_light))
            moon_terminator_line = termi[termi[:, 1] < cross_y, :]
            moon_terminator_line = np.vstack(([-cross_x,
                                               cross_y], moon_terminator_line))
            light = np.vstack((moon_edge_light, moon_terminator_line))
            moon_edge_shadow = edge[edge[:, 1] < cross_y, :]
            moon_edge_shadow = np.vstack(([cross_x,
                                           cross_y], moon_edge_shadow[::-1]))
            shadow = np.vstack((moon_terminator_line, moon_edge_shadow))
    else:
        light = edge[edge[:, 1] >= 0, :]
        shadow = edge[edge[:, 1] <= 0, :]

    # moon phase affine transformation
    M_rotation = np.array(
        [[np.cos(light_side_angle), -np.sin(light_side_angle)],
         [np.sin(light_side_angle),
          np.cos(light_side_angle)]])

    if light is not None:
        light = np.dot(M_rotation, light.T)
        light = moon_wcs.pixel_to_world(light[0], light[1])
        light = SkyCoord(
            light.ra, light.dec, frame=loc_gcrs
        )  # transform fake icrs coordinates back to location's gcrs
        light = light.transform_to(az_frame)  # convert to az coordinates
        light = np.array([light.az.value, 90 - light.alt.value]).T
        # convert altitude angle to zenith angle
        light = light / 180 * np.pi  # convert to rad
        light_codes = np.ones(len(light),
                              dtype=mpath.Path.code_type) * mpath.Path.LINETO
        light_codes[0] = mpath.Path.MOVETO
        light_path = mpath.Path(light, light_codes)
        light_patch = mpatches.PathPatch(light_path,
                                         facecolor=light_color,
                                         linewidth=0,
                                         zorder=999)
        ax.add_patch(light_patch)
    if shadow is not None:
        shadow = np.dot(M_rotation, shadow.T)
        shadow = moon_wcs.pixel_to_world(shadow[0], shadow[1])
        shadow = SkyCoord(shadow.ra, shadow.dec, frame=loc_gcrs)
        shadow = shadow.transform_to(az_frame)
        shadow = np.array([shadow.az.value, 90 - shadow.alt.value]).T
        shadow = shadow / 180 * np.pi
        shadow_codes = np.ones(len(shadow),
                               dtype=mpath.Path.code_type) * mpath.Path.LINETO
        shadow_codes[0] = mpath.Path.MOVETO
        shadow_path = mpath.Path(shadow, shadow_codes)
        shadow_patch = mpatches.PathPatch(shadow_path,
                                          facecolor=[shadow_color] * 3,
                                          linewidth=0,
                                          zorder=999)
        ax.add_patch(shadow_patch)

    # Lunar mare patch
    moon_now = get_moon(time=moon.obstime)
    moon_1h = get_moon(time=moon.obstime + 1 * u.h)
    # a very dirty way to calculate moon's north pole's postion angle, without considering moon's libration.
    pos_angle = position_angle(moon_now.ra, moon_now.dec, moon_1h.ra,
                                moon_1h.dec) - np.pi / 2 * u.rad
    # affine transformation for rotating by pos_angle
    NP_rotation = np.array([[np.cos(pos_angle), -np.sin(pos_angle)],
                            [np.sin(pos_angle),
                                np.cos(pos_angle)]])
    for lm in lms:
        lm = np.dot(NP_rotation, lm.T)
        lm = moon_wcs.pixel_to_world(lm[0], lm[1])
        lm = SkyCoord(lm.ra, lm.dec, frame=loc_gcrs)
        lm = lm.transform_to(az_frame)
        lm = np.array([lm.az.value, 90 - lm.alt.value]).T
        lm = lm / 180 * np.pi
        lm_codes = np.ones(len(lm),
                            dtype=mpath.Path.code_type) * mpath.Path.LINETO
        lm_codes[0] = mpath.Path.MOVETO
        lm_path = mpath.Path(lm, lm_codes)
        lm_patch = mpatches.PathPatch(lm_path,
                                        facecolor=[0.5, 0.5, 0.5, 0.6],
                                        linewidth=0,
                                        zorder=1000)
        ax.add_patch(lm_patch)


if __name__ == '__main__':
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz, GCRS, get_moon, get_body
    from astropy.time import Time
    import matplotlib.pyplot as plt
    import matplotlib.scale as mscale
    from stereographic import StereographicZenithScale
    mscale.register_scale(StereographicZenithScale)
    time = Time('2022-1-1 00:00:00')-6.742*u.h+2*u.day
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
    sun = get_body('sun', time, location)
    moon = get_moon(time, location)
    plt.figure()
    plt.subplot(111, projection="polar")
    plt.yscale('stereographiczenith')
    ax = plt.gca()
    ax.set_theta_zero_location("N")
    draw_moon_logo(ax, location, az, moon, sun, loc_gcrs, resolution=50)
    plt.show()