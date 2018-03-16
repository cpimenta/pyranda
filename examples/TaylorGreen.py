import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim
from pyrandaIBM import pyrandaIBM
from pyrandaBC  import pyrandaBC
from pyrandaTimestep  import pyrandaTimestep





problem = 'TGvortex'

## Define a mesh
Npts = 32
imesh = """
xdom = (0.0, 2*pi,  Npts, periodic=True)
ydom = (0.0, 2*pi,  Npts, periodic=True)
zdom = (0.0, 2*pi,  Npts, periodic=True)
""".replace('Npts',str(Npts)).replace('pi',str(numpy.pi))

    
# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,imesh)
ss.addPackage( pyrandaTimestep(ss) )


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)                  - ddz(:rho:*:w:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tauxx:) - ddy(:rhou:*:v: - :tauxy:)       - ddz(:rhou:*:v: - :tauxz:)
ddt(:rhov:) =  -ddx(:rhov:*:u: + :tauxy:)       - ddy(:rhov:*:v: + :p: - :tauyy:) - ddz(:rhov:*:w: - :tauyz:)
ddt(:rhow:) =  -ddx(:rhow:*:u: + :tauxz:)       - ddy(:rhow:*:v: + :tauyz:)       - ddz(:rhow:*:w: + :p: - :tauzz:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tauxx:)*:u: ) - ddy( (:Et: + :p: - :tauyy:)*:v: ) - ddz( (:Et: + :p: - :tauzz:)*:w: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity 
:ux:        =  ddx(:u:)
:vy:        =  ddy(:v:)
:wz:        =  ddz(:w:)
:div:       =  :ux: + :vy: + :wz:
# Remaining cross derivatives
:uy:        =  ddy(:u:)
:uz:        =  ddz(:u:)
:vx:        =  ddx(:v:)
:vz:        =  ddz(:v:)
:wy:        =  ddy(:w:)
:wx:        =  ddx(:w:)
:enst:      = sqrt( (:uy:-:vx:)**2 + (:uz: - :wx:)**2 + (:vz:-:wy:)**2 )
:tke:       = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + :wz:*:wz: + .5*((:uy:+:vx:)**2 + (:uz: + :wx:)**2 + (:vz:+:wy:)**2) )
:mu:        =  gbar( ring(:S:  ) ) * :rho: * 1.0e-4
:beta:      =  gbar( ring(:div:) ) * :rho: * 0.0e-2
:tauxx:     =  2.0*:mu:*:ux:   + (:beta:-2./3.*:mu:) *:div:
:tauyy:     =  2.0*:mu:*:vy:   + (:beta:-2./3.*:mu:) *:div:
:tauzz:     =  2.0*:mu:*:wz:   + (:beta:-2./3.*:mu:) *:div:
:tauxy:     = :mu:*(:uy:+:vx:) + (:beta:-2./3.*:mu:) *:div:
:tauxz:     = :mu:*(:uz:+:wx:) + (:beta:-2./3.*:mu:) *:div:
:tauyz:     = :mu:*(:vz:+:wz:) + (:beta:-2./3.*:mu:) *:div:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
#:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
"""

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
u0 = 1.0
p0 = 1.0
rho0 = 1.0
L = 1.0
:u: =  u0*sin(meshx/L)*cos(meshy/L)*cos(meshz/L)
:v: = -u0*cos(meshx/L)*sin(meshy/L)*cos(meshz/L)
:w: = 0.0*:u:
:p:  = p0 + rho0/16.0*( cos(2.*meshx/L) + cos(2.*meshy/L))*(cos(2.*meshz/L) + 2.0)
:rho: = rho0 + 0.0*:u:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:Et:  = :p: / (:gamma:-1.0) + 0.5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:tke: = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

# Set the initial conditions
ss.setIC(ic)
    
# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]


# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time


# Mesh for viz on master
xx   =  ss.PyMPI.zbar( x ) / Npts
yy   =  ss.PyMPI.zbar( y ) / Npts
#xx   =   x[:,:,16] / Npts
#yy   =   y[:,:,16] / Npts
ny = ss.PyMPI.ny

# Start time loop
CFL = 0.5
dt = ss.variables['dt'].data * CFL

# Viz
cnt = 1
viz_freq = 50
pvar = 'u'

tke0 = ss.variables['tke'].data.sum()

while time < 20.0:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = ss.variables['dt'].data * CFL

    # Print some output
    maxU = ss.variables['u'].data.max()
    tke = ss.variables['tke'].data.sum()/tke0
    #ss.iprint("%s -- %s --- Max u: %s" % (cnt,time,maxU)  )
    ss.iprint("%s -- %s --- TKE: %s" % (cnt,time,tke)  )
    cnt += 1
    if viz:
        v = ss.PyMPI.zbar( ss.variables[pvar].data )
        #v = ss.variables[pvar].data[:,:,16]
        if (ss.PyMPI.master and (cnt%viz_freq == 0)) and True:
            plt.figure(2)
            plt.clf()            
            plt.contourf( xx,yy,v ,64 , cmap=cm.jet)
            plt.title(pvar)
            plt.pause(.001)



