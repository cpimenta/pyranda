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
        pdat = self.pyranda.variables[var].pydata
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                for dd in range(pdat.rank):
                    pdat.data[dd][0,:,:] = pdat.data[dd][1,:,:]

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                for dd in range(pdat.rank):
                    pdat.data[dd][-1,:,:] = pdat.data[dd][-2,:,:]

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                for dd in range(pdat.rank):
                    pdat.data[dd][:,0,:] = pdat.data[dd][:,1,:]

        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                for dd in range(pdat.rank):
                    pdat.data[dd][:,-1,:] = pdat.data[dd][:,-2,:]

                
    def const(self,var,direction,val):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.constant( v , d , val)
                
                
    def constant(self,var,direction,val):

        pdat = self.pyranda.variables[var].pydata
        
        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                for dd in range(pdat.rank):
                    pdat.data[dd][0,:,:] = val
                
            
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                for dd in range(pdat.rank):
                    pdat.data[dd][-1,:,:] = val
            

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                for dd in range(pdat.rank):
                    pdat.data[dd][:,0,:] = val
                
            
        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                for dd in range(pdat.rank):
                    pdat.data[dd][:,-1,:] = val
            
            
        


