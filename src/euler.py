from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')
import matplotlib.pyplot as plt
from pyranda import pyrandaSim,pyrandaMPI,fortran3d

from ibm import pyrandaIBM


## Define a mesh
Npts = 128
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts

mesh_options = {}
mesh_options['type'] = 'cartesian'
mesh_options['periodic'] = numpy.array([False, True, True])
mesh_options['dim'] = 3
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp   , Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, 1 ,  1  ]
#mesh_options['nn'] = [ Npts, Npts ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',mesh_options)
ss.addPackage( pyrandaIBM(ss) )


# Define the equations of motion
eom ="""
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)   - ddy(:rhou:*:v:)
ddt(:rhov:) =  -ddx(:rhov:*:u:)                 - ddy(:rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: ) - ddy( (:Et: + :p: - :tau:)*:v: )
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(lap(lap(:div:))))*:dx: * :rho: * 0.2
:tau:       =  :beta:*:div:
[:gx:,:gy:,:gz:] = grad( :phi: )
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:] )
:Et:          = ibmS( :Et:  , :phi:, [:gx:,:gy:,:gz:] )
:rho:         = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( 1.4 - 1.0 )
"""
ss.EOM(eom)

#:umag:      =  numpy.sqrt(:u:*:u: + :v:*:v:)

# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]


rad = numpy.sqrt( (x-numpy.pi)**2 ) # + (y-numpy.pi)**2 ) #+ (z-numpy.pi)**2  )
ss.variables['Et'].data = 1.0 + 0.1 * numpy.exp( -(rad)**2/(.2**2) )
ss.variables['rho'].data += 1.0
#ss.variables['p'].data += 1.0 
ss.variables['dx'].data += (x[1,0,0] - x[0,0,0])**6

ss.variables['Et'].data =  ss.gfilter(numpy.where( x < numpy.pi , 1.0 / .4, .1 / .4))
#ss.variables['Et'].data =  ss.gfilter(numpy.where( x < numpy.pi , 1.0 / .4, 0.99999 / .4))
ss.variables['rho'].data = ss.gfilter(numpy.where( x < numpy.pi , 1.0, .125))
#ss.variables['rho'].data = ss.gfilter(numpy.where( x < numpy.pi , 1.0, 1.0))
ss.variables['dx'].data += (x[1,0,0] - x[0,0,0])**6

ss.variables['phi'].data = 4.0 - x


time = 0.0
viz = True

#import pdb
#pdb.set_trace()

v = 1.0

dt_max = v / ss.mesh.nn[0] * 1.0

tt = L/v * .5 #dt_max

xx   =  ss.PyMPI.zbar( x )
yy   =  ss.PyMPI.zbar( y )

dt = dt_max
cnt = 1


while tt > time:

    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:
        phi = ss.PyMPI.zbar( ss.variables['u'].data )
        if ss.PyMPI.master and (cnt%40 == 0) and True:
            plt.figure(2)
            plt.clf()
            #plt.contourf( xx,yy,phi ,64 )
            plt.plot(xx[:,0],phi[:,0] ,'b.-')
            plt.pause(.001)



