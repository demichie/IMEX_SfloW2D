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

if len(sys.argv)==3: 
 
    runname = sys.argv[1]

    anim_duration = int(sys.argv[2])

else:

    print('Please provide two arguments:\n')
    print('1) Run name (as given in IMEX_SfloW2D.inp)\n')
    print('2) Anim duration in seconds\n')
    sys.exit()

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

plt.xlim([np.amin(x_unique),np.amax(x_unique)])

line1, = ax.plot(x_unique,B_cent)
line2, = ax.plot(x_unique,hB)  

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


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
