#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
from mpl_toolkits.mplot3d import Axes3D                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys
import os.path
import matplotlib.tri as mtri
from matplotlib import cm


if len(sys.argv)==2: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 

else:

    print('Please provide two or three arguments:\n')
    print('1) File name\n')
    print('2) 1st Variable to plot: h,hB,B,u,v\n')
    print('2) 2nd Variable to plot: h,hB,B,u,v\n')
    sys.exit()

data = np.loadtxt(filename,skiprows=0)

x = data[:,0]
y = data[:,1]
h = data[:,2]
u = data[:,3]
v = data[:,4]
b = data[:,5]
w = data[:,6]

x0_idx = np.asarray((np.where(data[:,0]==data[0,0])))

ny_cells = x0_idx[0,1]
nx_cells = data.shape[0] / ny_cells

X_cent = x.reshape((nx_cells,ny_cells))
Y_cent = y.reshape((nx_cells,ny_cells))
H_cent = h.reshape((nx_cells,ny_cells))
B_cent = b.reshape((nx_cells,ny_cells))
W_cent = w.reshape((nx_cells,ny_cells))

W_nan = np.copy(W_cent)
W_nan[H_cent<1e-10] = np.nan 

# create a figure for the plot
fig = plt.figure()
ax = fig.gca(projection='3d')


surf = ax.plot_surface(X_cent, Y_cent, W_nan,linewidth=0, antialiased=False, vmin=-1, vmax=1)
ax.plot_surface(X_cent, Y_cent, W_cent, alpha=0.3)
cset = ax.contour(X_cent, Y_cent, H_cent, np.logspace(-4, -1, num=4),zdir='z', linewidths=(0.5,),offset=1000, cmap=cm.coolwarm)

ax.set_xlabel('X')
#ax.set_xlim(0, 30)
ax.set_ylabel('Y')
#ax.set_ylim(-15, 15)
ax.set_zlabel('Z')
ax.set_zlim(1000, 3500)
#ax.axis('equal')
ax.view_init(11,-33)
# fig.tight_layout()
plt.show()    


