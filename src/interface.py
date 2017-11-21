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
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u: - :tau:)                 
ddt(:rhoA:)  =  -ddx(:rhoA:*:uA: - :tauA:)  
ddt(:phi:)  =  - :gx: * :u1:
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhoA:      =  fbar( :rhoA:  )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho: * 0.0 + 0.7
:uA:        =  :rhoAuA: / :rhoA: * 0.0 + 0.7
:dr:        = ddx(:rho:)
:drA:       = ddx(:rhoA:)
:tau:       = 4.0e-5*gbar(abs(lap(lap(:rho:)))) * :dx4: * :u: * :dr: * :dx:
:tauA:      = 4.0e-5*gbar(abs(lap(lap(:rhoA:)))) * :dx4: * :uA: * :drA: * :dx:
# Immersed boundary method used for MM levelset
[:gx:,:gy:,:gz:] = grad( :phi: )
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,0.0,0.0] )
:rho:         = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
[:uA:,:v:,:w:] = ibmV( [:uA:,:v:,0.0], -:phi:, [-:gx:,-:gy:,-:gz:], [:u1:,0.0,0.0] )
:rhoA:         = ibmS( :rhoA: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:uT:        =  where( :phi: > 0.0, :u:, :uA:)
:rhoB:        =  where( :phi: > 0.0, :rho:, :rhoA:)
:rhoT:        =  where( abs(:phi:/:dx:*2.0) > 1.0, :rhoB:, .5*(:phi:/:dx:*2.0 + 1.0)*(:rho:-:rhoA:) + (:rhoA:)   )
"""
ss.EOM(eom)

#:umag:      =  numpy.sqrt(:u:*:u: + :v:*:v:)

# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]


rad = numpy.sqrt( (x-numpy.pi)**2 ) # + (y-numpy.pi)**2 ) #+ (z-numpy.pi)**2  )
#ss.variables['rho'].data += 1.0 + 1.0e6* (1.0-(1.0 + numpy.tanh( (x-numpy.pi*.25)/0.08 ))*.5)
ss.variables['rho'].data += ss.gfilter( numpy.where( x < numpy.pi*.25, 1.0e6, 1.0 ) )

ss.variables['rhoA'].data += 1.0e6 * (1.0 + .1*numpy.exp(-(x-numpy.pi*1.5)**2/.01)  )

# Init momenta
ss.variables['u1'].data = .70

ss.variables['dx'].data += (x[1,0,0] - x[0,0,0])
ss.variables['dx4'].data += (x[1,0,0] - x[0,0,0])**4

# Define the interface location
ss.variables['phi'].data = 3.14159 - x


time = 0.0
viz = True

v = 1.0

dt_max = v / ss.mesh.nn[0] * .15

tt = L/v * .25 #dt_max

xx   =  ss.PyMPI.zbar( x )
yy   =  ss.PyMPI.zbar( y )

dt = dt_max
cnt = 1


while tt > time:

    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )
    ss.variables['rho'].data[0,:,:] = 1.0e6 + 1.0
    
    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:
        #v1 = ss.PyMPI.zbar( ss.variables['rhoA'].data )
        v2 = ss.PyMPI.zbar( ss.variables['rhoT'].data )
        if ss.PyMPI.master and (cnt%5 == 0) and True:
            plt.figure(2)
            plt.clf()
            #plt.contourf( xx,yy,phi ,64 )
            #plt.plot(xx[:,0],v1[:,0] ,'b.-')
            plt.plot(xx[:,0],v2[:,0] ,'k.-')
            plt.pause(.001)
            #import pdb
            #pdb.set_trace()



