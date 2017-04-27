"""
maxdist

Calculate distance limits in midplane for various surveys.
Request of Blanton for SDSS-IV overview paper.

"""
from __future__ import print_function
import numpy as np
from scipy.interpolate import UnivariateSpline
import os
import pandas as pd
os.environ['DUST_DIR'] = '/Users/jquark/obs_data/dustmaps/'
import mwdust

Hdustmap = mwdust.Combined15(filter='2MASS H')
Vdustmap = mwdust.Combined15(filter='CTIO V')
distances = np.linspace(.1, 50, 100)

solar_2gyr_isochrone = pd.read_csv('parsec_iso.txt', sep="\t", header=13)
mhpoor_6_2gyr_isochrone = pd.read_csv('parsec_z0.00352', sep="\t", header=13)
# Get 1.765 Msol RC star ; minus 1 from index for lower mass
solar_ind = np.searchsorted(solar_2gyr_isochrone['M_ini'], 1.765) - 1
solar_source_rc = solar_2gyr_isochrone.loc[int(solar_ind)]
solar_source_rg = solar_2gyr_isochrone.loc[125]
solar_source_tiprgb = solar_2gyr_isochrone.loc[136]

metalpoor_source_rg = mhpoor_6_2gyr_isochrone.loc[140]

dustmaps = {}
dustmaps['H'] = Hdustmap
dustmaps['V'] = Vdustmap


def dist2distmod(dist):
    """dist in kpc to distance modulus"""
    return 5.*np.log10(dist)+10.


def distmod2dist(distmod):
    """distance modulus to distance in kpc"""
    return 10.**(distmod/5.-2.)


def los_maxdist(l, b, band='H', source_mag=-3.0, mag_limit=11.0):
    """
    Given array of distances along a LOS in Galaxy, calculate max distance
    observable given maximum distance modulus (depends on survey and
    absolute magnitude of source.

    Parameters
    ----------------

    distances : array
        set of increasing distances along LOS
    Av : array
        set of extinctions corresponding to distances
    maxDM : float
        maximum distance modulus, m - M, where m is the limiting magnitude
        of your survey and M is the absolute magnitude of your source in
        the same photometric band as your extinction.

    Returns
    ----------------
    distance

    """
    # absMag = rc_mag[band]
    distmod = mag_limit - source_mag
    Av = dustmaps[band](l, b, distances)
    print(l, b, band)
    RHS = ((distmod - 10. - Av)/5.)
    LHS = np.log10(distances)
    fit = UnivariateSpline(distances, RHS-LHS, s=0)
    return fit.roots()[0]  # distance in kiloparsecs

apogee = {'band': 'H', 'mag_limit': 12.2}
galah = {'band': 'V', 'mag_limit': 14.0}
gaia_eso = {'band': 'V', 'mag_limit': 19.0}


def survey_dist(survey, source_mag=-5.0):
    glon = np.linspace(0, 360.0, 100.0)
    glat = np.zeros_like(glon)
    maxdist = np.zeros_like(glon)
    for ii, (l, b) in enumerate(zip(glon, glat)):
        maxdist[ii] = los_maxdist(l, b, band=survey['band'],
                                  source_mag=source_mag,
                                  mag_limit=survey['mag_limit'])

    results = {}
    results['l'] = glon
    results['b'] = glat
    results['maxdist'] = maxdist
    return results


def save_maxdist(filename, source_star=solar_source_rg):
    apodata = survey_dist(apogee, source_mag=source_star.H)
    galahdata = survey_dist(galah, source_mag=source_star.V)
    gaia_esodata = survey_dist(gaia_eso, source_mag=source_star.V)
    results = {}
    results['l'] = apodata['l']
    results['b'] = apodata['b']
    results['apogee'] = apodata['maxdist']
    results['galah'] = galahdata['maxdist']
    results['gaia_eso'] = gaia_esodata['maxdist']
    df = pd.DataFrame.from_dict(results)
    df.to_pickle(filename+'.df')
    df.to_csv(filename+'.txt', sep='\t', columns=['l', 'b', 'apogee', 'galah', 'gaia_eso'])

if __name__ == '__main__':
    save_maxdist('limits_solarRC', source_star=solar_source_rc)
    save_maxdist('limits_solarRG', source_star=solar_source_rg)
    save_maxdist('limits_solarTipRGB', source_star=solar_source_tiprgb)
