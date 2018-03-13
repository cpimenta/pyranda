import sys
import time
import numpy 
import matplotlib.pyplot as plt
from pyranda import pyrandaSim


# Define the domain/mesh
domain = "xdom = (0.0 , 1.0 , 100 )"

# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',domain)

# Define the equations of motion
ss.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")

# Initialize variables
ic = """
:phi: = 1.0 + 0.1 * exp( -(abs(:x:-.5)/.1 )**2 )
:phi0: = :phi:
:c:   = 1.0
"""
ss.setIC(ic)

# Integrate in time
dt = .001 
time = 0.0
while time < .1:
    time = ss.rk4(time,dt)


# Plot the initial/final solution
x   = ss.mesh.coords[0]
phi = ss.variables['phi'].data
phi0 = ss.variables['phi0'].data
plt.plot(x[:,0,0],phi[:,0,0] ,'k.-')
plt.plot(x[:,0,0],phi0[:,0,0] ,'b-')
plt.show()




            


