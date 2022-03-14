from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS, GCRS, get_moon, get_body
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.scale as mscale
from stereographic import StereographicZenithScale

mscale.register_scale(StereographicZenithScale)

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
    The ellipse of moon's terminator. return an array of its pixel coordinates [X,Y] with 2 column. The points are not uniform in the arc.
    '''
    t = np.linspace(0, np.pi * 2.0, resolution)
    t = t.reshape((len(t), 1))
    x = np.cos(t) * k
    y = np.sin(t) / b * k
    return np.hstack((x, y))


def draw_moon_logo(az_frame,
                   moon,
                   sun,
                   loc_gcrs,
                   zoom_factor=20,
                   resolution=50,
                   light_color=0.75,
                   shadow_color=0.25,
                   k=1.005):
    '''
    Get the matplotlib patch of moon's light and shadow in the latitude-zenith coordinate. Return in tuple
    '''
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
    cross_x = np.sqrt((b**2 - k**2) / (b**2 - 1))
    cross_y = np.sqrt((k**2 - 1) / (b**2 - 1))
    if np.cos(moon_phase) > 0:
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

    M_rotation = np.array(
        [[np.cos(light_side_angle), -np.sin(light_side_angle)],
         [np.sin(light_side_angle),
          np.cos(light_side_angle)]])
    light = np.dot(M_rotation, light.T)
    shadow = np.dot(M_rotation, shadow.T)

    light = moon_wcs.pixel_to_world(light[0], light[1])
    light = SkyCoord(
        light.ra, light.dec, frame=loc_gcrs
    )  # transform fake icrs coordinates back to location's gcrs
    light = light.transform_to(az_frame)  # convert to az coordinates
    light = np.array(
        [light.az.value, 90 - light.alt.value]
    ).T  # convert back to the path format at the same time convert altitude angle to zenith angle
    light = light / 180 * np.pi  # convert to rad

    shadow = moon_wcs.pixel_to_world(shadow[0], shadow[1])
    shadow = SkyCoord(shadow.ra, shadow.dec, frame=loc_gcrs)
    shadow = shadow.transform_to(az_frame)
    shadow = np.array([shadow.az.value, 90 - shadow.alt.value]).T
    shadow = shadow / 180 * np.pi

    light_codes = np.ones(len(light),
                          dtype=mpath.Path.code_type) * mpath.Path.LINETO
    light_codes[0] = mpath.Path.MOVETO
    light_path = mpath.Path(light, light_codes)
    light_patch = mpatches.PathPatch(light_path,
                                     facecolor=[light_color] * 3,
                                     linewidth=0,
                                     zorder=999)

    shadow_codes = np.ones(len(shadow),
                           dtype=mpath.Path.code_type) * mpath.Path.LINETO
    shadow_codes[0] = mpath.Path.MOVETO
    shadow_path = mpath.Path(shadow, shadow_codes)
    shadow_patch = mpatches.PathPatch(shadow_path,
                                      facecolor=[shadow_color] * 3,
                                      linewidth=0,
                                      zorder=999)

    return (light_patch, shadow_patch)


if __name__ == '__main__':
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

    sun = get_body('sun', time, location)
    moon = get_moon(time, location)

    r = np.linspace(30, 60, 50) / 180 * np.pi
    theta = np.linspace(0, 2 * np.pi, 50)
    plt.figure()
    plt.subplot(111, projection="polar")
    plt.plot(
        theta,
        r,
        '-',
        lw=2,
    )
    plt.yscale('stereographiczenith')
    ax = plt.gca()
    ax.set_theta_zero_location("N")
    l_p, s_p = draw_moon_logo(az, moon, sun,resolution = 50)
    ax.add_patch(l_p)
    ax.add_patch(s_p)
    plt.show()