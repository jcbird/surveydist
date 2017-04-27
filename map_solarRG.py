from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from galpy.util import bovy_coords


solar_2gyr_isochrone = pd.read_csv('parsec_iso.txt', sep="\t", header=13)
# Get 1.765 Msol RC star ; minus 1 from index for lower mass
solar_source_rg = solar_2gyr_isochrone.loc[125]

#  Coordinates for Galactic Disk (radius= 15 kpc)
galaxy_rad = 15.0
x_gc = np.linspace(-galaxy_rad, galaxy_rad, 100)
y_gc = np.sqrt(galaxy_rad**2 - x_gc**2)
minusy_gc = -1.0*y_gc
x_gc_solarpos = x_gc + 8.0  # GC is at x=8.0

#  Concatenate for full circle
galaxy_xs = np.concatenate((x_gc_solarpos[::-1], x_gc_solarpos))
galaxy_ys = np.concatenate((y_gc[::-1], minusy_gc))
galaxy_zs = np.zeros_like(galaxy_xs)
galaxy_lbd = bovy_coords.XYZ_to_lbd(galaxy_xs, galaxy_ys, galaxy_zs)

#  Coordinates for Galactic Bulge (radius= 3 kpc)
galaxy_rad = 3.0
x_bulge = np.linspace(-galaxy_rad, galaxy_rad, 100)
y_bulge = np.sqrt(galaxy_rad**2 - x_bulge**2)
minusy_bulge = -1.0*y_bulge
x_bulge_solar = x_bulge + 8.0  # GC is at x=8.0
#  Concatenate for full circle
bulge_xs = np.concatenate((x_bulge_solar[::-1], x_bulge_solar))
bulge_ys = np.concatenate((y_bulge[::-1], minusy_bulge))
bulge_zs = np.zeros_like(bulge_xs)
bulge_lbd = bovy_coords.XYZ_to_lbd(bulge_xs, bulge_ys, bulge_zs)

#
bulge_l = np.linspace(-np.arctan(3.5/8.), np.arctan(3.5/8.), 100)
coeffs = [np.repeat(1, len(bulge_l)), -2*8.0*np.cos(bulge_l),
          np.repeat(8**2-3.0**2, len(bulge_l))]
bulge_ds = np.array([np.roots(ii) for ii in np.array(coeffs).T])

labels = [r'APOGEE (m$_{\mathrm{H}}=12.2$)', r'GALAH (m$_{\mathrm{V}}=14.0$)',
          r'GAIA-ESO (m$_{\mathrm{V}}=19.0$)']


def plot_maxdist(maxdist_df_file, filename='test'):
    max_distances = pd.read_pickle(maxdist_df_file)
    with plt.style.context(['ggplot', 'avrfig']):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='polar')
        for ii, survey in enumerate(['apogee', 'galah', 'gaia_eso']):
            baseline = ax.plot(np.radians(max_distances['l']),
                               max_distances[survey].values, label=labels[ii])
            fill = ax.fill_between(np.radians(max_distances['l']),
                                   np.zeros_like(max_distances['l']),
                                   max_distances[survey].values,
                                   facecolor=baseline[0].get_color(),
                                   alpha=.15)
        # ax.bar(np.radians(345.0), 11.0, width=np.radians(30.0), bottom=5.0)
        # ax.plot(galaxy_lbd[:,0], galaxy_lbd[:,2], lw=1, c='k')
        bulge = ax.plot(bulge_lbd[:, 0],  bulge_lbd[:, 2],  lw=1,
                        label='MW Bulge', c='k')
        ax.fill_between(bulge_l, bulge_ds[:, 0], bulge_ds[:, 1],
                        facecolor=bulge[0].get_color(), alpha=.15)
        ax.set_rlabel_position(180.0)
        ax.set_rmax(15.0)
        plt.legend(bbox_to_anchor=(1.3, 1.05), fontsize=13, loc=1)
        #  filename = 'DistLim_1.72Msun_solarZ_RG'
        plt.savefig(filename + '.png', format='png')
        plt.savefig(filename + '.pdf', format='pdf')
        # ax.set_rasterized(True)
        # plt.savefig(filename + '.eps', format='eps')
        # ax.plot(galaxy_lbd[:,0], galaxy_lbd[:,2], lw=1, c='k')
        # ax.plot(galaxy_lbd[:,0], galaxy_lbd[:,2], lw=1, c='k')

if __name__ == '__main__':
    plot_maxdist('limits_solarRC.df', filename='solarRC_map')
    plot_maxdist('limits_solarRG.df', filename='solarRG_map')
    plot_maxdist('limits_solarTipRGB.df', filename='solarTipRGB_map')
