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
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(lap(lap(:div:))))*:dx6: * :rho: * 0.2
:tau:       =  :beta:*:div:
"""


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = "rad = sqrt( (:x:-numpy.pi)**2  ) "
if dim == 2:
    ic = "rad = sqrt( (:x:-numpy.pi)**2  +  (:y:-numpy.pi)**2 ) "

# Linear wave propagation in 1d and 2d
if (problem == 'linear'):
    pvar = 'p'
    ic += """
    :gamma: = 1.4
    ratio = 1.0 + 0.01 * exp( -(rad)**2/(.2**2) )
    :Et: = ratio
    :rho: = 1.0
    """

# SOD shock tube in 1d and 2d
if (problem == 'sod'):
    pvar = 'rho'
    if dim == 1:
        ic = 'rad = :x: / 2.0'
    ic += """
    :gamma: = 1.4
    :Et:  = gbar( where( rad < :pi:/2.0, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
    :rho: = gbar( where( rad < :pi:/2.0, 1.0    , .125 ) )
    """

# Set the initial conditions
ss.setIC(ic)
    
# Length scale for art. viscosity
x = ss.mesh.coords[0]
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



