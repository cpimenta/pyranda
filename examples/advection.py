#from mpi4py import MPI
import numpy 
import re
import sys
import time
#sys.path.append('/Users/olson45/Research/FloATPy')
sys.path.append('../')
sys.path.append('../../python_tools/compact-light')
sys.path.append('../../python_tools/compac-light/mpi4py/install/lib/python2.7/site-packages')


import matplotlib.pyplot as plt
from pyranda.pyranda import pyrandaSim

from pyranda.pyrandaIBM import pyrandaIBM


# Try to get args
try:
    Npts = int(sys.argv[1])
except:
    Npts = 100

try:
    test = bool(sys.argv[2])
except:
    test = False

## Define a mesh
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts


mesh_options = {}
mesh_options['type'] = 'cartesian'
mesh_options['periodic'] = numpy.array([True, True, True])
mesh_options['dim'] = 1
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp   , Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, 1 ,  1  ]

# Initialize a simulation object on a mesh
#import pdb
#pdb.set_trace()
ss = pyrandaSim('advection',mesh_options)

# Define the equations of motion
ss.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")

# Initialize variables
ic = """
rad   = sqrt( (:x:-:pi:)**2  )
:phi: = 1.0 + 0.1 * exp( -(rad/(:pi:/4.0))**2 )
:phi2: = :phi:*1.0
:c:   = 1.0
"""
ss.setIC(ic)

x  = ss.mesh.coords[0]
xx =  ss.PyMPI.zbar( x )

# Time step size
v = 1.0
dt_max = v / ss.mesh.nn[0] * L * .90
tt = L/v * 1.0 


# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = True
while tt > time:

    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    if not test:
        ss.iprint("%s -- %s" % (cnt,time)  )

    # Plot animation of advection
    cnt += 1
    if viz:
        v = ss.PyMPI.zbar( ss.variables['phi'].data )
        if (ss.PyMPI.master and (cnt%5 == 0)) and (not test):
            plt.figure(1)
            plt.clf()
            plt.plot(xx[:,0],v[:,0] ,'k.-')
            plt.pause(.001)


phi = ss.variables['phi'].data
phi2 = ss.variables['phi2'].data
error = numpy.sum( (phi-phi2)**2  )
ss.iprint( error ) 
            


