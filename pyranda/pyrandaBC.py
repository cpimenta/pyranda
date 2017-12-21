import numpy 
import re
import sys
import time
from pyrandaPackage import pyrandaPackage


class pyrandaBC(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'BC'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.bcList = {}
        

    def get_sMap(self):
        sMap = {}
        sMap['bc.extrap('] = "self.packages['BC'].extrap("
        sMap['bc.const('] =  "self.packages['BC'].const("
        self.sMap = sMap
        

    
        
    def extrap(self,var,direction):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.extrapolate( v , d )
            
        

    def extrapolate(self,var,direction):
        # Direction switch
        bcvar = None
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = self.pyranda.variables[var].data[1,:,:]

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = self.pyranda.variables[var].data[-2,:,:]

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = self.pyranda.variables[var].data[:,1,:]

        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = self.pyranda.variables[var].data[:,-2,:]

                
    def const(self,var,direction,val):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.constant( v , d , val)
                
                
    def constant(self,var,direction,val):

        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = val
                
            
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = val
            

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = val
                
            
        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = val
            
            
        


