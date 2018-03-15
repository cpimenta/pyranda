from mpi4py import MPI
import numpy 
import re
import sys
import time


#sys.path.append('/Users/olson45/Research/FloATPy')
#from floatpy.parallel import t3dmod
#from floatpy.derivatives.compact import CompactDerivative
#from floatpy.filters.filter import Filter                
            
import parcop


class pyrandaMPI():

    #def __init__(self,nx,ny,nz,dx,dy,dz,periodic):
    def __init__(self,meshOptions):

        self.nx = meshOptions['nn'][0]
        self.ny = meshOptions['nn'][1]
        self.nz = meshOptions['nn'][2]
        
        x1 = meshOptions['x1'][0]
        xn = meshOptions['xn'][0]
        y1 = meshOptions['x1'][1]
        yn = meshOptions['xn'][1]
        z1 = meshOptions['x1'][2]
        zn = meshOptions['xn'][2]
        
        dx = (meshOptions['xn'][0]-meshOptions['x1'][0])/max(self.nx-1,1)
        dy = (meshOptions['xn'][1]-meshOptions['x1'][1])/max(self.ny-1,1)
        dz = (meshOptions['xn'][2]-meshOptions['x1'][2])/max(self.nz-1,1)
        self.dx = dx
        self.dy = dy
        self.dz = dz
        
        periodic = meshOptions['periodic']
        self.periodic = periodic
        
        self.comm  = MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
            
            
        self.order = (10,10,10)
        self.filter_type = ('compact', 'compact', 'compact')
        
        px = 1
        py = 1
        pz = 1
        bx1 = "NONE"
        bxn = "NONE"
        by1 = "NONE"#"PERI"
        byn = "NONE"#"PERI"
        bz1 = "NONE"#"PERI"
        bzn = "NONE"#"PERI"
        if periodic[0]:
            bx1 = "PERI"
            bxn = "PERI"
        if periodic[1]:
            by1 = "PERI"
            byn = "PERI"
        if periodic[2]:
            bz1 = "PERI"
            bzn = "PERI"
        

        parcop.parcop.setup( 0, 0 , self.fcomm, self.nx,self.ny,self.nz, 
                             px , py, pz, x1,xn,y1,yn,
                             z1,zn,bx1,bxn,by1,byn,
                             bz1,bzn)
        
            #self.grid_partition = t3dmod.t3d(self.fcomm,
            #                                 self.nx,self.ny,self.nz,
            #                                 self.periodic)

        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')
        
            #self.grid_partition.get_sz3d(self.chunk_3d_size)
            #self.grid_partition.get_st3d(self.chunk_3d_lo)
            #self.grid_partition.get_en3d(self.chunk_3d_hi)

            
        # Set up comms
        #self.xcom  = MPI.Comm.f2py( self.grid_partition.commx () )
        #self.ycom  = MPI.Comm.f2py( self.grid_partition.commy () )
        #self.zcom  = MPI.Comm.f2py( self.grid_partition.commz () )
        #self.xycom = MPI.Comm.f2py( self.grid_partition.commxy() )
        #self.yzcom = MPI.Comm.f2py( self.grid_partition.commyz() )
        #self.xzcom = MPI.Comm.f2py( self.grid_partition.commxz() )
        self.xcom  = MPI.Comm.f2py( parcop.parcop.commx (0,0) ) 
        self.ycom  = MPI.Comm.f2py( parcop.parcop.commy (0,0) )
        self.zcom  = MPI.Comm.f2py( parcop.parcop.commz (0,0) )
        self.xycom = MPI.Comm.f2py( parcop.parcop.commxy(0,0) )
        self.yzcom = MPI.Comm.f2py( parcop.parcop.commyz(0,0) )
        self.xzcom = MPI.Comm.f2py( parcop.parcop.commxz(0,0) )

            

        self.chunk_3d_size[0] = self.nx / px
        self.chunk_3d_size[1] = self.ny / py
        self.chunk_3d_size[2] = self.nz / pz

        self.chunk_3d_lo[0] = self.xcom.rank     * self.chunk_3d_size[0] + 1
        self.chunk_3d_hi[0] = (self.xcom.rank+1) * self.chunk_3d_size[0]

        self.chunk_3d_lo[1] = self.ycom.rank     * self.chunk_3d_size[1] + 1
        self.chunk_3d_hi[1] = (self.ycom.rank+1) * self.chunk_3d_size[1]
            
        self.chunk_3d_lo[2] = self.zcom.rank     * self.chunk_3d_size[2] + 1
        self.chunk_3d_hi[2] = (self.zcom.rank+1) * self.chunk_3d_size[2]

        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing



            # Replace these objects with parcop ones
            #self.der = CompactDerivative(self.grid_partition,
            #(self.dx, self.dy, self.dz),
            #                             self.order, self.periodic)
            #self.fil = Filter( self.grid_partition, self.filter_type, periodic_dimensions=self.periodic )
            #self.gfil = Filter( self.grid_partition, ('gaussian','gaussian','gaussian'), periodic_dimensions=self.periodic ) 

        self.der  = parcop_der()
        self.fil  = parcop_sfil()
        self.gfil = parcop_gfil()


        self.patch = 0
        self.level = 0
        self.setPatch()        

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
                                 
    def setPatch(self):                                 
        parcop.parcop.set_patch( self.patch, self.level )

    def emptyScalar(self):
        return numpy.zeros( self.chunk_3d_size, dtype=numpy.float64, order='F')*0.0

    def emptyVector(self):
        blk_size = numpy.append(self.chunk_3d_size,3)
        return numpy.zeros( blk_size, dtype=numpy.float64, order='F')*0.0

    
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
    




class parcop_der:

    def __init__(self):        
        pass

    def ddx(self,val):        
        return parcop.parcop.ddx( val )

    def ddy(self,val):        
        return parcop.parcop.ddy( val )

    def ddz(self,val):        
        return parcop.parcop.ddz(  val )

    def laplacian(self,val):        
        return parcop.parcop.plaplacian(  val )

    def ring(self,val):
        return parcop.parcop.pring(  val )

class parcop_gfil:

    def __init__(self):
        pass

    def filter(self,val):
        return parcop.parcop.gfilter( val)

class parcop_sfil:

    def __init__(self):
        pass

    def filter(self,val):
        return parcop.parcop.sfilter(  val)
