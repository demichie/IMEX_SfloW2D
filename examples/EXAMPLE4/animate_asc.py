#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
from linecache import getline
import sys
import os.path
import time

from mayavi import mlab
import  moviepy.editor as mpy

if len(sys.argv)==3: 
 
    runname = sys.argv[1]
    anim_length = sys.argv[2]

else:

    print('Please provide two arguments:\n')
    print('1) run name (as given in IMEX_SfloW2D.inp)\n')
    print('2) duration of the animated gif in seconds\n')
    sys.exit()

bakfile = runname+'.bak'

with open(bakfile) as f:
    for line in f:
        if "DT_OUTPUT" in line:
             line1 = line.replace(',','')
             dt = float(line1.split()[-1])
             print 'Sim dt='+str(dt)
        if "T_END" in line:
             line1 = line.replace(',','')
             duration = float(line1.split()[-1])
             print 'Sim duration='+str(duration)


fps_orig = 1.0/dt
anim_length = 5.0
speed_factor = duration / anim_length 


source1 = 'dem_esri.asc'
source2 = runname+'_{0:04}'.format(0)+'.asc'

# Parse the header using a loop and
# the built-in linecache module
hdr = [getline(source1, i) for i in range(1,7)]
values = [float(h.split(" ")[-1].strip()) \
 for h in hdr]
cols,rows,lx,ly,cell,nd = values
xres = cell
yres = cell * -1

# Load the dem into a numpy array
arr = np.loadtxt(source1, skiprows=6)

nx = arr.shape[1]
xs = lx +0.5*cell + np.linspace(0,(nx-1)*cell,nx)
xmin = np.min(xs)
xmax = np.max(xs)

ny = arr.shape[0]
ys = ly+0.5*cell - np.linspace(0,(ny-1)*cell,ny)
ymin = np.min(ys)
ymax = np.max(ys)

ys = np.linspace(ymin,ymax,ny)

B_cent = np.zeros((ny,nx))

for i in range(0,ny):

   B_cent[i,0:nx] = arr[i,0:nx]

H_cent = np.zeros_like(B_cent)
W_cent = np.zeros_like(B_cent)

# MAKE A FIGURE WITH MAYAVI

fig_myv = mlab.figure(size=(600,600),bgcolor=(1,1,1))
# fig_myv = mlab.figure(bgcolor=(1,1,1))

# ANIMATE THE FIGURE WITH MOVIEPY, WRITE AN ANIMATED GIF

def make_frame(t):
    i = int(np.ceil(t*speed_factor/dt))
    source2 = runname+'_{0:04}'.format(i)+'.asc'
    # Load the dem into a numpy array
    arr = np.loadtxt(source2, skiprows=6)

    for i in range(0,ny):

        H_cent[i,0:nx] = arr[i,0:nx]

    
    W_cent = B_cent + H_cent
    idx = np.ma.masked_where(H_cent==-9999,W_cent)
    W_cent[np.where(np.ma.getmask(idx)==True)] = np.nan

   
    mlab.clf() # clear the figure (to reset the colors)
    topo = mlab.surf(xs,ys, B_cent,color=(0.5,0.6,0.7))
    mlab.outline(topo, color=(.7,.7,.7))
    surf = mlab.surf(xs,ys, W_cent+1.e-1,color=(0.4,0.4,0.4))
    surf.actor.property.interpolation = 'phong'
    surf.actor.property.specular = 1.0
    surf.actor.property.specular_power = 50
    mlab.text(-10,0,'t='+str(t*speed_factor)+'s',z=10,width=0.15,color=(0,0,0))
    return mlab.screenshot(antialiased=True)


animation = mpy.VideoClip(make_frame, duration=duration/speed_factor)
animation.write_gif(runname+".gif", fps=fps_orig*speed_factor)
animation.write_videofile(runname+".mp4", fps=fps_orig*speed_factor)




