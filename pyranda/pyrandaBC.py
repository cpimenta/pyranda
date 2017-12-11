import numpy 
import re
import sys
import time
from pyrandaPackage import pyrandaPackage


class pyrandaBC(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'BC'
        pyrandaPackage.__init__(self,PackageName,pysim)


    def get_sMap(self):
        sMap = {}
        sMap['bc.extrap('] = "self.packages['BC'].extrapolate("
        self.sMap = sMap
        

    def extrapolate(self,var,direction):

        # Direction switch
        bcvar = None
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                bcvar = self.pyranda.variables[var].data[1,:,:]

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                bcvar = self.pyranda.variables[var].data[-2,:,:]


        return bcvar


                

        
    def BCconstant(var,direction,val):

        # Direction switch
        if direction == 'x1':
            if self.PyMPI.xcom.rank == 0:
                self.variables[var].data[0,:,:] = val
                
            
        if direction == 'xn':
            if self.PyMPI.xcom.rank == self.PyMPI.xcom.size - 1:
                self.variables[var].data[-1,:,:] = val
            
            

            
        


