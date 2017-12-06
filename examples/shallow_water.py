from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')
import matplotlib.pyplot as plt
from pyranda import pyrandaSim,pyrandaMPI


## Define a mesh
Npts = 64
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts

mesh_options = {}
mesh_options['type'] = 'cartesian'
mesh_options['periodic'] = numpy.array([True, True, True])
mesh_options['dim'] = 3
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp   , Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, Npts ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',mesh_options)


# Define the equations of motion
eom ="""
ddt(:eta:)  = - ddx(:eta:*:u:) - ddy(:eta:*:v:)
ddt(:ueta:) = - ddx(:ueta:*:u: + .5*:g:*:eta:**2) - ddy(:ueta:*:v:)
ddt(:veta:) = - ddx(:veta:*:u:) - ddy(:veta:*:v: + .5*:g:*:eta:**2)
:eta:       =  fbar( :eta:  )
:ueta:      =  fbar( :ueta: )
:veta:      =  fbar( :veta: )
:v:         =  :veta: / :eta:
:u:         =  :ueta: / :eta:
"""
ss.EOM(eom)

# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]

rad = numpy.sqrt( (x-numpy.pi)**2 + (y-numpy.pi)**2 ) #+ (z-numpy.pi)**2  )
ss.variables['eta'].data = 1.0 + .01 * numpy.exp( -(rad)**2/(.2**2) )
ss.variables['g'].data = 1.0


time = 0.0
viz = True


v = 1.0

dt_max = v / ss.mesh.nn[0] * .75

tt = L/v * .5 #dt_max

phi1 = ss.PyMPI.zbar( ss.variables['eta'].data )
xx   =  ss.PyMPI.zbar( x )
yy   =  ss.PyMPI.zbar( y )
if ss.PyMPI.master:
    plt.figure(1)
    plt.contour( xx,yy,phi1 , 32 )
    plt.figure(2)
    plt.plot(xx[:,Npts/2],phi1[:,Npts/2],'k--')
    


dt = dt_max
cnt = 1
while tt > time:

    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:
        phi = ss.PyMPI.zbar( ss.variables['eta'].data )
        if ss.PyMPI.master and (cnt%10 == 0):
            plt.figure(2)
            plt.clf()
            plt.contourf( xx,yy,phi ,32 )
            plt.pause(.001)


phi = ss.PyMPI.zbar( ss.variables['eta'].data )
if ss.PyMPI.master:
    plt.figure(1)
    plt.contour( xx,yy,phi , 32 )
    plt.figure(2)
    plt.plot(xx[:,Npts/2],phi[:,Npts/2],'b-')

    plt.figure(3)
    plt.plot(xx[:,Npts/2],phi[:,Npts/2] ,'b-')

    plt.show()
