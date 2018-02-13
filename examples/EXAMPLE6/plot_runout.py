#!/usr/bin/env python
"""
% This function 

"""
import numpy as np
from linecache import getline
from mpl_toolkits.mplot3d import Axes3D                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import matplotlib.tri as mtri

import time
import sys

filename = sys.argv[1]

# create a figure for the plot
fig, ax = plt.subplots()

t,dist = np.loadtxt(filename,usecols=(3,7), unpack=True)

plt.xlim([np.amin(t),np.amax(t)])

line1, = ax.plot(t,dist)

ax.set_xlabel('time (s)')
ax.set_ylabel('runout (m)')

plt.show()
