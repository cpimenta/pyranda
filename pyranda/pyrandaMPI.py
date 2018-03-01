from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')

from floatpy.parallel import t3dmod
from floatpy.derivatives.compact import CompactDerivative
from floatpy.filters.filter import Filter                
            
class pyrandaMPI():

    def __init__(self,nx,ny,nz,dx,dy,dz,periodic):

            self.nx = nx
            self.ny = ny
            self.nz = nz

            self.dx = dx
            self.dy = dy
            self.dz = dz

            self.comm  = MPI.COMM_WORLD
            self.fcomm = self.comm.py2f()
            self.periodic = periodic

            self.order = (10,10,10)
            self.filter_type = ('compact', 'compact', 'compact')

            self.grid_partition = t3dmod.t3d(self.fcomm,
                                             self.nx,self.ny,self.nz,
                                             self.periodic)

            self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
            self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
            self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
            self.grid_partition.get_sz3d(self.chunk_3d_size)
            self.grid_partition.get_st3d(self.chunk_3d_lo)
            self.grid_partition.get_en3d(self.chunk_3d_hi)
            self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
            self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing

            # Set up comms
            self.xcom  = MPI.Comm.f2py( self.grid_partition.commx () )
            self.ycom  = MPI.Comm.f2py( self.grid_partition.commy () )
            self.zcom  = MPI.Comm.f2py( self.grid_partition.commz () )
            self.xycom = MPI.Comm.f2py( self.grid_partition.commxy() )
            self.yzcom = MPI.Comm.f2py( self.grid_partition.commyz() )
            self.xzcom = MPI.Comm.f2py( self.grid_partition.commxz() )
            
            self.der = CompactDerivative(self.grid_partition,
                                         (self.dx, self.dy, self.dz),
                                         self.order, self.periodic)

            self.fil = Filter( self.grid_partition, self.filter_type, periodic_dimensions=self.periodic )
            self.gfil = Filter( self.grid_partition, ('gaussian','gaussian','gaussian'), periodic_dimensions=self.periodic ) 
            self.master = False
            if self.comm.rank == 0:
                self.master = True


            self.x1proc = False
            if self.xcom.rank == 0:
                self.x1proc = True

            self.xnproc = False
            if self.xcom.rank == self.xcom.size - 1:
                self.xnproc = True

            self.y1proc = False
            if self.ycom.rank == 0:
                self.y1proc = True

            self.ynproc = False
            if self.ycom.rank == self.ycom.size - 1:
                self.ynproc = True


    # Data initializers for scalars, vectors, and tensors
    def emptyScalar(self):
        return numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F')*0.0

    def emptyVector(self):
        return [ self.emptyScalar(), self.emptyScalar(), self.emptyScalar() ]

    def emptyTensor(self):
        return [ self.emptyVector(), self.emptyVector(), self.emptyVector() ]
    
    

    
    def gather2D(self,idata,com,n1,n2,g1,g2):
        """
        Give idata(a1,a2) return the assembled
          global array of gdata(n1,n2)
        """        

        ldata = numpy.zeros( (n1,n2) )
        ldata[ g1[0]:g1[1], g2[0]:g2[1] ] = idata

        fldata = ldata.flatten()
        fldataT = fldata * 0.0
        com.Allreduce( fldata, fldataT, op=MPI.SUM )
        
        return numpy.reshape(fldataT, (n1,n2), order = 'C' )


    def sum2D(self,data,com,n2,n3,index,g2,g3):

        a2 = g2[1] - g2[0]
        a3 = g3[1] - g3[0]
        gsum = numpy.zeros((a2,a3))

        # Get the local proc mean
        lsum = numpy.sum( data, index )

        tsum = self.gather2D(lsum,com,n2,n3,g2,g3)

        return tsum

    def sum3D(self,data):

        lsum = numpy.sum( data )
        tsum = self.comm.allreduce( lsum, op=MPI.SUM )

    def max3D(self,data):

        lmax = numpy.max( data )
        tmax = self.comm.allreduce( lmax, op=MPI.MAX )

        return tmax

    def min3D(self,data):

        lmin = numpy.max( data )
        tmin = self.comm.allreduce( lmin, op=MPI.MIN )

        return tmin

    def xbar(self,data):        
        return self.sum2D( data, self.comm,self.ny,self.nz,0,
                           [ self.chunk_3d_lo[1] , self.chunk_3d_hi[1] ],
                           [ self.chunk_3d_lo[2] , self.chunk_3d_hi[2] ]  )
  
    def ybar(self,data):
        return self.sum2D( data, self.comm,self.nx,self.nz,1,
                           [self.chunk_3d_lo[0] , self.chunk_3d_hi[0]],
                           [self.chunk_3d_lo[2] , self.chunk_3d_hi[2]])
    
    def zbar(self,data):

        return self.sum2D( data, self.comm,self.nx,self.ny,2,
                           [self.chunk_3d_lo[0],self.chunk_3d_hi[0]+1 ],
                           [self.chunk_3d_lo[1],self.chunk_3d_hi[1]+1 ])
    
