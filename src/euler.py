from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')
import matplotlib.pyplot as plt
from matplotlib import cm
from pyranda import pyrandaSim,pyrandaMPI,fortran3d

from ibm import pyrandaIBM


## Define a mesh
Npts = 128
L = numpy.pi * 2.0  
dim = 2
gamma = 1.4

problem = 'linear'
problem = 'sod'

Lp = L * (Npts-1.0) / Npts
mesh_options = {}
mesh_options['type'] = 'cartesian'
mesh_options['periodic'] = numpy.array([False, True, True])
mesh_options['dim'] = 3
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp   , Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, 1 ,  1  ]
if dim == 2:
    mesh_options['nn'] = [ Npts, Npts ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',mesh_options)


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)   - ddy(:rhou:*:v:)
ddt(:rhov:) =  -ddx(:rhov:*:u:)                 - ddy(:rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: ) - ddy( (:Et: + :p: - :tau:)*:v: )
ddt(:phi:)  =  - :gx: * :u1:
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( %s - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(lap(lap(:div:))))*:dx6: * :rho: * 0.2
:tau:       =  :beta:*:div:
""" % ( gamma )


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]

rad = numpy.sqrt( (x-numpy.pi)**2  )
if dim == 2:
    rad = numpy.sqrt( (x-numpy.pi)**2  + (y-numpy.pi)**2 )
    
if (problem == 'linear'):
    pvar = 'p'
    ratio = 1.0 + 0.01 * numpy.exp( -(rad)**2/(.2**2) )
    ss.variables['Et'].data = 1.0*ratio
    ss.variables['rho'].data += 1.0 
    
if (problem == 'sod'):
    pvar = 'rho'
    if dim == 1:
        rad = x/2.0
    ss.variables['Et'].data =  ss.gfilter(numpy.where(rad < numpy.pi/2.0 , 1.0 / .4,.1 / .4))
    ss.variables['rho'].data = ss.gfilter(numpy.where(rad < numpy.pi/2.0 , 1.0, .125))

# Length scale for art. viscosity
ss.variables['dx6'].data += (x[1,0,0] - x[0,0,0])**6


# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.5
tt = L/v * .25 #dt_max

# Mesh for viz on master
xx   =  ss.PyMPI.zbar( x )
yy   =  ss.PyMPI.zbar( y )

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 25
ss.updateVars()
while tt > time:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    
    # Print some output
    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:
        v = ss.PyMPI.zbar( ss.variables[pvar].data )
        if (ss.PyMPI.master and (cnt%viz_freq == 1)) and True:
            plt.figure(1)
            plt.clf()
            ny = ss.PyMPI.ny
            if ( ny > 1):
                plt.plot(xx[:,ny/2],v[:,ny/2] ,'k.-')
                plt.title(pvar)
                plt.pause(.001)
                plt.figure(2)
                plt.clf()            
                plt.contourf( xx,yy,v ,64 , cmap=cm.jet)
            else:
                plt.plot(xx[:,0],v[:,0] ,'k.-')
            plt.title(pvar)
            plt.pause(.001)



