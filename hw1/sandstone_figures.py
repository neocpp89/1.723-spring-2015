#!/usr/bin/env python
import matplotlib
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rc_file(r'../mpl.rc')
matplotlib.rc('lines', linewidth=2)

with open('berea_xsection_bot.dat', 'r') as fbot:
    bot_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), fbot)
    bot_data = np.array(bot_data)

    plt.figure(figsize=(4,4))
    plt.title('Berea Sandstone (bottom surface)')
    #plt.xlabel('X [pixels]')
    #plt.ylabel('Y [pixels]')
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BW', matplotlib.cm.get_cmap('Greys',10)(np.arange(10)))
    norm = matplotlib.colors.BoundaryNorm([0,0.5,1], cmap.N)
    plt.imshow(bot_data, cmap=cmap, norm=norm)
    plt.axis('off')
    cb = plt.colorbar()
    cb.ax.set_yticklabels(['0', '', '1'])
    plt.gca().set_aspect('equal')
    plt.savefig('bot.pdf')

with open('berea_xsection_top.dat', 'r') as ftop:
    top_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), ftop)
    top_data = np.array(top_data)

    plt.figure(figsize=(4,4))
    plt.title('Berea Sandstone (top surface)')
    #plt.xlabel('X [pixels]')
    #plt.ylabel('Y [pixels]')
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('BW', matplotlib.cm.get_cmap('Greys',10)(np.arange(10)))
    norm = matplotlib.colors.BoundaryNorm([0,0.5,1], cmap.N)
    plt.imshow(top_data, cmap=cmap, norm=norm)
    plt.axis('off')
    cb = plt.colorbar()
    cb.ax.set_yticklabels(['0', '', '1'])
    plt.gca().set_aspect('equal')
    plt.savefig('top.pdf')
