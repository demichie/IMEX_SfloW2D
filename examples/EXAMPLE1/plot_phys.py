#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys
import os.path

print 

if len(sys.argv)==3: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 

    var = sys.argv[2]

elif len(sys.argv)==4: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 

    var = sys.argv[2]

    var2 = sys.argv[3]

else:

    print('Please provide two or three arguments:\n')
    print('1) File name\n')
    print('2) 1st Variable to plot: h,hB,B,u,v\n')
    print('2) 2nd Variable to plot: h,hB,B,u,v\n')
    sys.exit()

data = np.loadtxt(filename,skiprows=0)

print data.shape

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

print n_cent,n_unique

# create a figure for the plot
fig, ax = plt.subplots()


if ( n_cent == n_unique):

    if var=='h':

        z = h

    elif var=='B':

        z = B_cent

    elif var=='u':

        z = u

    elif var=='hu':

        z = h*u

    elif var=='v':

        z = v

    elif var=='hB':

        z = hB

    else:

        print('Please specify the variable to plot as 2nd argument: h,hB,B,u,v')
        sys.exit()


    plt.xlim([np.amin(x_unique),np.amax(x_unique)])

    line1, = ax.plot(x_unique,z)

    if len(sys.argv)==4:

        if var2=='h':

            z2 = h

        elif var2=='B':

            z2 = B_cent

        elif var2=='u':

            z2 = u

        elif var2=='v':

            z2 = v

        elif var2=='hB':

            z2 = hB

        else:

            print('Please specify the 2nd variable to plot as 3nd argument: h,hB,B,u,v')
            sys.exit()

        line2, = ax.plot(x_unique,z2)




plt.show()





