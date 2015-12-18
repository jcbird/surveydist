import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import os

os.environ['DUST_DIR']='/hd1/obs_data/dustmaps/'


thetas = np.radians(np.linspace(0, 359.9, 40))
rads = np.random.random(40)+2.0

with plt.style.context(['ggplot', 'avrfig']):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    fill_apogee = ax.fill_between(thetas, np.zeros_like(rads), rads, zorder=0,
                              alpha=.15, label='APOGEE')
    apogee_color = fill_apogee.get_facecolor()[0]
    ax.plot([], [], color=apogee_color[:3], alpha=apogee_color[-1], 
            label='APOGEE', lw=10)
    ax.set_rlabel_position(210.0)
    plt.legend(bbox_to_anchor=(1.05, 1.05), loc=1)
    plt.show()
