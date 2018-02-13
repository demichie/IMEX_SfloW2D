#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys

if len(sys.argv)==4: 

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        n_cells = int(a)
        print n_cells
 
    else:
 
        sys.exit()

else:

    print('Please provide three arguments:\n')
    print('1) Number of cells\n')
    print('2) Variables to reconstruct: phys or cons\n')
    print('3) Order of the RK scheme\n')
    sys.exit()

# Define the boundaries x_left and x_right of the spatial domain
x_left = -0.25
x_right = 1.75
y_bottom = 0.0

# Define the number n_points of points of the grid
n_points  = n_cells+1

# Define the grid stepsize dx
dx = ( x_right - x_left ) / ( n_points - 1 )


# Define the array x of the grid points
x = np.linspace(x_left,x_right,n_points)

x_cent = np.linspace(x_left+0.5*dx,x_right-0.5*dx,n_cells)

y = np.zeros((n_points,1))

y_cent = np.zeros_like(x_cent)
w_cent = np.zeros_like(x_cent)
u_cent = np.zeros_like(x_cent)



# define the topography
for i in range(n_points):

    if x[i]<0:
    
        y[i] = 1
    
    elif x[i]<=0.4:
    
        y[i] = (np.cos(np.pi*x[i]))**2
    
    elif x[i]<=0.5:
    
        y[i] = (np.cos(np.pi*x[i]))**2 + 0.25 * ( np.cos(10.0*np.pi*(x[i]-0.5))+1)
    
    elif x[i]<= 0.6:
    
        y[i] = 0.5*(np.cos(np.pi*x[i]))**4 + 0.25 * ( np.cos(10.0*np.pi*(x[i]-0.5))+1)

    elif x[i]<= 1.0:
    
        y[i] = 0.5*(np.cos(np.pi*x[i]))**4
    
    elif x[i]<= 1.5:
    
        y[i] = 0.25*(np.sin(2.0*np.pi*(x[i]-1.0)))
    
    else:
    
        y[i] = 0.0
    

# define the initial solution
for i in range(n_cells):

    y_cent[i] = 0.5*(y[i]+y[i+1])

    if x_cent[i]<0:
    
        w_cent[i] = 1.4
        u_cent[i] = 0.0
    
    else: 
 
        w_cent[i] = y_cent[i]
        u_cent[i] = 0.0



# create topography file
header = "ncols     %s\n" % n_points
header += "nrows    %s\n" % 1
header += "xllcorner " + str(x_left-0.5*dx) +"\n"
header += "yllcorner " + str(0-0.5*dx) +"\n"
header += "cellsize " + str(dx) +"\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, np.transpose(y), header=header, fmt='%1.12f',comments='')

# create intial solution file
q0 = np.zeros((5,n_cells))

q0[0,:] = x_cent
q0[1,:] = 0.0
q0[2,:] = w_cent
q0[3,:] = (w_cent-y_cent)*u_cent
q0[4,:] = 0.0

np.savetxt('example3_0000.q_2d', np.transpose(q0), fmt='%19.12e') 

with open('example3_0000.q_2d','a') as file:
    file.write('\n')

# Read in the file
with open('IMEX_SfloW2D.template', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'example3')
filedata = filedata.replace('restartfile', 'example3_0000.q_2d')
filedata = filedata.replace('x_left', str(x_left))
filedata = filedata.replace('y_bottom', str(y_bottom))
filedata = filedata.replace('n_cells', str(n_cells))
filedata = filedata.replace('dx', str(dx))

if sys.argv[2]=='cons':

    print('Linear reconstruction of conservative variables (h+B,hu,hv)')
    filedata = filedata.replace('bc2', 'HU')
    filedata = filedata.replace('recvar', 'cons')

else:

    print('Linear reconstruction of physical variables (h+B,u,v)')
    filedata = filedata.replace('bc2', 'U')
    filedata = filedata.replace('recvar', 'phys')

filedata = filedata.replace('order', sys.argv[3])

# Write the file out again
with open('IMEX_SfloW2D.inp', 'w') as file:
  file.write(filedata)

# create a figure for the plot
fig, ax = plt.subplots()
plt.ylim([-0.1,1.5])
plt.xlim([-0.25,1.75])

# plot the initial solution and call "line" the plot
line1, = ax.plot(x,y)
line2, = ax.plot(x_cent,w_cent,'-o')


plt.show()
