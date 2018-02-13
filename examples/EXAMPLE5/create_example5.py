#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys

def isfloat(string):
  try:
    return float(string) and '.' in string  # True if string is a number contains a dot
  except ValueError:  # String is not a number
    return False


if len(sys.argv)==6: 

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        n_cells = int(a)
        print n_cells
 
    else:
 
        sys.exit()

    print('Topography Slope (in degrees)')
    a = sys.argv[2]

    if isfloat(a) or a.isdigit():

        slope = float(a)
        print slope
 
    else:
 
        sys.exit()

    print('Pile Slope')
    a = sys.argv[3]

    if isfloat(a) or a.isdigit():

        pileslope = float(a)
        print pileslope
 
    else:
 
        sys.exit()


else:

    print('Please provide three arguments:\n')
    print('1) Number of cells\n')
    print('2) Topography Slope (in degrees)\n') 
    print('3) Pile Slope\n')
    print('4) Variables to reconstruct: phys or cons\n')
    print('5) Order of the RK scheme\n')
    sys.exit()


# Define the boundaries x_left and x_right of the spatial domain
x_left = 0.0
x_right = 500.0


# Define the number n_points of points of the grid
n_points  = n_cells+1

# Define the grid stepsize dx
dx = ( x_right - x_left ) / ( n_points - 1 )


# Define the array x of the grid points
x = np.linspace(x_left,x_right,n_points)

x_cent = np.linspace(x_left+0.5*dx,x_right-0.5*dx,n_cells)

B = np.zeros_like(x)

B_cent = np.zeros_like(x_cent)
w_cent = np.zeros_like(x_cent)
u_cent = np.zeros_like(x_cent)



# define the topography
B[0:n_points] = ( x_right - x[0:n_points] ) * np.tan( slope/180.0 * np.pi )  

#B[0:n_points] = 0.003*( x[0:n_points] -250.0 ) **2  


# define the initial solution
for i in range(n_cells):

    B_cent[i] = 0.5*(B[i]+B[i+1])

    if (x_cent[i]>200) and ( x_cent[i] < 300 ):
    
        w_cent[i] = B_cent[i] + (50.0-np.abs(x_cent[i]-250)) * np.tan( pileslope/180.0 * np.pi )  
        #w_cent[i] = B_cent[i] + 20.0
        u_cent[i] = 0.0
    
    else: 
 
        w_cent[i] = B_cent[i]
        u_cent[i] = 0.0



# create topography file
header = "ncols     %s\n" % n_points
header += "nrows    %s\n" % 1
header += "xllcorner " + str(x_left-0.5*dx) +"\n"
header += "yllcorner " + str(0-0.5*dx) +"\n"
header += "cellsize " + str(dx) +"\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, np.transpose(B), header=header, fmt='%1.12f',comments='')

# create intial solution file
q0 = np.zeros((6,n_cells))

q0[0,:] = x_cent
q0[1,:] = 0.0
q0[2,:] = w_cent
q0[3,:] = (w_cent-B_cent)*u_cent
q0[4,:] = 0.0
q0[5,:] = (w_cent-B_cent)*0.2

np.savetxt('example5_0000.q_2d', np.transpose(q0), fmt='%19.12e') 

with open('example5_0000.q_2d','a') as file:
    file.write('\n')

# Read in the file
with open('IMEX_SfloW2D.template', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'example5')
filedata = filedata.replace('restartfile', 'example5_0000.q_2d')
filedata = filedata.replace('x_left', str(x_left))
filedata = filedata.replace('x_right', str(0.0))
filedata = filedata.replace('n_cells', str(n_cells))
filedata = filedata.replace('dx', str(dx))

if sys.argv[4]=='cons':

    print('Linear reconstruction of conservative variables (h+B,hu,hv)')
    filedata = filedata.replace('bc2', 'HU')
    filedata = filedata.replace('recvar', 'cons')

else:

    print('Linear reconstruction of physical variables (h+B,u,v)')
    filedata = filedata.replace('bc2', 'U')
    filedata = filedata.replace('recvar', 'phys')

filedata = filedata.replace('order', sys.argv[5])


# Write the file out again
with open('IMEX_SfloW2D.inp', 'w') as file:
  file.write(filedata)

# create a figure for the plot
fig, ax = plt.subplots()
fig1 = plt.gcf()


dx = x_cent[-1] - x_cent[0]
dz = np.amax(w_cent) - np.amin(B_cent)

fig_ratio = dz/dx

fig1.set_size_inches(10/fig_ratio, 10, forward=True)
plt.axis('equal')
plt.xlim([0,500])

# plot the initial solution and call "line" the plot
line1, = ax.plot(x,B)
line2, = ax.plot(x_cent,w_cent,'-')


plt.show()
