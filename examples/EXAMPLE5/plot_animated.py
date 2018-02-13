#!/usr/bin/env python

"""
Animate the 1D output of IMEX-SfloW2D

"""

import glob
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys
import os.path

runname = sys.argv[1]

anim_duration = int(sys.argv[2])

bakfile = runname+'.bak'

with open(bakfile) as f:
    for line in f:
        if "DT_OUTPUT" in line:
             line1 = line.replace(',','')
             dt = float(line1.split()[-1])
             print 'dt='+str(dt)
        if "T_END" in line:
             line1 = line.replace(',','')
             duration = float(line1.split()[-1])
             print 'Sim duration='+str(duration)



filelist = sorted(glob.glob(runname+'_*[0-9].p_2d'))

filename = runname+'_{0:04}'.format(0)+'.p_2d'
data = np.loadtxt(filename,skiprows=0)

x_cent = data[:,0]
y_cent = data[:,1]
h = data[:,2]
u = data[:,3]
v = data[:,4]
B_cent = data[:,5]
hB = data[:,6]

x_unique = np.unique(x_cent)

n_cent = len(x_cent)
n_unique = len(x_unique)

# create a figure for the plot
fig, ax = plt.subplots()
fig1 = plt.gcf()

dx = x_cent[-1] - x_cent[0]
dz = np.amax(hB) - np.amin(B_cent)

fig_ratio = dz/dx

#fig1.set_size_inches(10/fig_ratio, 10, forward=True)
plt.axis('equal')
plt.xlim([0,500])

line1, = ax.plot(x_unique,B_cent)
line2, = ax.plot(x_unique,hB)  

time_template = 'time = %.1fs'
time_text = ax.text(0.80, 0.9, '', transform=ax.transAxes)


ax.grid()


def animate(i):

    filename = runname+'_{0:04}'.format(i)+'.p_2d'
    # print(filename)
    data = np.loadtxt(filename,skiprows=0)

    x_cent = data[:,0]
    hB = data[:,6]

    x_unique = np.unique(x_cent)


    line2.set_data(x_unique,hB)  
    time_text.set_text(time_template % (i*dt))
    
    return line2,time_text

n_frames = len(filelist)
print 'n_frames=',n_frames

millisec = anim_duration * 1000.0 / ( n_frames - 1) 
print 'millisec=',millisec

anim_fps = ( n_frames - 1 )/ anim_duration
print 'anim_fps=',anim_fps

ani = animation.FuncAnimation(fig, animate, np.arange(0, len(filelist)),
                              interval=millisec)

ani.save(runname+'.mp4', fps=anim_fps)
ani.save(runname+'.gif', fps=anim_fps, writer='imagemagick')
plt.show()
