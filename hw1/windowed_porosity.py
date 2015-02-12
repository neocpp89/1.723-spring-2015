#!/usr/bin/env python
import matplotlib
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import sys
matplotlib.rc_file(r'../mpl.rc')
matplotlib.rc('lines', linewidth=2)

from calculate_porosity import calculate_porosity

def calculate_expanding_porosity_window(data):
    num_rows = len(data)
    num_cols = len(data[0, :])
    # assume square
    if (num_cols != num_rows):
        sys.exit(0)
    window_edge_size = []
    porosities = []
    for i in range(0, num_rows/2):
        s = i
        lens = num_rows - 2*i
        window_edge_size.append(lens)
        t = calculate_porosity(data[s:(s+lens),s:(s+lens)])
        porosities.append(t[2])

    window_edge_size.reverse()
    porosities.reverse()

    #print window_edge_size
    #print porosities

    return (window_edge_size, porosities)

with open('berea_xsection_bot.dat', 'r') as fbot:
    data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), fbot)
    data = np.array(data)
    bot_windowed_porosity = calculate_expanding_porosity_window(data)

with open('berea_xsection_top.dat', 'r') as ftop:
    data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), ftop)
    data = np.array(data)
    top_windowed_porosity = calculate_expanding_porosity_window(data)

plt.figure(figsize=(4,4))
plt.title('Porosity versus Window Size')
plt.xlabel(r'Window Edge Size [$\sqrt{\#\mathrm{\ pixels}}$]')
plt.ylabel('Porosity` [-]')
plt.plot(bot_windowed_porosity[0], bot_windowed_porosity[1], ':k')
plt.plot(top_windowed_porosity[0], top_windowed_porosity[1], '-b')
#plt.ylim([0, 1])
plt.legend({'Bottom Surface','Top Surface'})
plt.tight_layout(pad=0.05)
plt.savefig('windowed_porosity.pdf')
