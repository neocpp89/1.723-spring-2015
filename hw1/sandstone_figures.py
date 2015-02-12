#!/usr/bin/env python
import matplotlib
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc_file(r'../mpl.rc')
matplotlib.rc('lines', linewidth=2)

def plot_sandstone_image(data, title, savename):
    plt.figure(figsize=(3.5,3.5))
    plt.title(title)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BW', matplotlib.cm.get_cmap('Greys',10)(np.arange(10)))
    norm = matplotlib.colors.BoundaryNorm([0,0.5,1], cmap.N)
    plt.imshow(data, cmap=cmap, norm=norm)
    plt.axis('off')
    cb = plt.colorbar()
    cb.ax.set_yticklabels(['Void', '', 'Solid'])
    plt.gca().set_aspect('equal')
    plt.savefig(savename)
    return

if __name__ == "__main__":
    with open('berea_xsection_bot.dat', 'r') as fbot:
        data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), fbot)
        data = np.array(data)
        plot_sandstone_image(data, 'Berea Sandstone (bottom surface)', 'bot.pdf')

    with open('berea_xsection_top.dat', 'r') as ftop:
        data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), ftop)
        data = np.array(data)
        plot_sandstone_image(data, 'Berea Sandstone (top surface)', 'top.pdf')
