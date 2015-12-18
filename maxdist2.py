"""
maxdist

Calculate distance limits in midplane for various surveys.
Request of Blanton for SDSS-IV overview paper.

"""
from __future__ import print_function
import pandas as pd
import numpy as np
import mwdust
# from dustmap_query import query


_RvAv_file = 'extinction_table.txt'


def ex_from_reddening(band=None):
    """
    Summary

    Parameters
    ----------------

    band : string
            Name of photometric band to calculate extinction from
            reddening. Uses extinction_table.txt

    Returns
    ----------------

    """
    conv_table = pd.read_csv('extinction_table.txt', sep='\t', comment='#',
                             names=['band', 'lambda', 'rv2.1', 'rv3.1', 'rv4.1', 'rv5.1'],
                             usecols=[0, 1, 2, 3, 4, 5])
    conv_table = conv_table.set_index('band')
    try:
        return conv_table.loc[band].loc['rv3.1']
    except ValueError:
        print('Function argument must be one of the following')
        print(conv_table.index.values)

def dist2distmod(dist):
    """dist in kpc to distance modulus"""
    return 5.*numpy.log10(dist)+10.

def distmod2dist(distmod):
    """distance modulus to distance in kpc"""
    return 10.**(distmod/5.-2.)

def calc_maxdist(l=30.0, b=0.0, source_Mag=2.0, limiting_mag=17.0, band=None):
    mapresult = query(l, b, coordsys='gal')
    distmod = limiting_mag - source_Mag
    #  mapresult distances are in distance moduli
    dm_index = np.searchsorted(mapresult['distmod'], distmod)
    if mapresult['converged'] is not 1:
        print('Caution, Not converged!')
    EBminusV = mapresult['best'][dm_index]
    conv_EBV_A = ex_from_reddening(band=band)
    extinction = EBminusV * conv_EBV_A
    maxdist = 10.**((distmod - 10.0 - extinction)/5.)
    return maxdist








