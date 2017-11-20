from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')

from floatpy.parallel import t3dmod
from floatpy.derivatives.compact import CompactDerivative
from floatpy.filters.filter import Filter

import matplotlib.pyplot as plt

USE_LAMBDA = True

class pyrandaMesh:

    def __init__(self,mesh_options,pympi):

        self.name = 'base'
        self.kind = mesh_options['type']
        self.options = mesh_options
        self.dims = mesh_options['dim']
        self.PyMPI = pympi
        self.x = None
        self.y = None
        self.z = None
        self.coords = None
        self.shape = None
        self.makeMesh()
        

    def makeMesh(self):


        options = self.options
        if self.kind == 'cartesian':

            x1 = options['x1']
            xn = options['xn']
            nn = options['nn']

            self.nn = nn    

            chunk_lo = self.PyMPI.chunk_3d_lo
            chunk_hi = self.PyMPI.chunk_3d_hi
            x = numpy.linspace(x1[0],xn[0],num=nn[0]+1)[chunk_lo[0]:chunk_hi[0]+1]
            y = numpy.linspace(x1[1],xn[1],num=nn[1]+1)[chunk_lo[1]:chunk_hi[1]+1]
            z = numpy.linspace(x1[2],xn[2],num=nn[2]+1)[chunk_lo[2]:chunk_hi[2]+1]

            x = numpy.linspace(x1[0],xn[0],num=nn[0])[chunk_lo[0]:chunk_hi[0]+1]
            y = numpy.linspace(x1[1],xn[1],num=nn[1])[chunk_lo[1]:chunk_hi[1]+1]
            z = numpy.linspace(x1[2],xn[2],num=nn[2])[chunk_lo[2]:chunk_hi[2]+1]
            
            x, y, z = numpy.meshgrid(x, y, z, indexing='ij')

            self.coords = [0]*3
            self.coords[0] = numpy.asfortranarray( x )
            self.coords[1] = numpy.asfortranarray( y )
            self.coords[2] = numpy.asfortranarray( z )

            self.shape = list(self.coords[0].shape)

            
            

class pyrandaVar:

    def __init__(self,name,kind,rank):

        self.name = name
        self.kind = kind
        self.rank = rank

        self.rhs = None   # RHS for algebraic closure
        self.dt  = None   # dt for temporal integration
        
    def allocate(self,pympi):
        
        # Inherit the mesh size
        if self.rank == 'scalar':
            self.data = pympi.emptyScalar()
        else:
            raise ValueError('Error: rank: %s not recognized' % self.rank)
            
        
class pyrandaEq:

    def __init__(self,eqstr,sMap):
        """
        Read in eqstr and extract the LHS variable 
        and the kind of equation (PDE or ALG)
        """
        self.eqstr = eqstr
        self.kind = 'ALG'
        self.LHS = findVar(eqstr.split('=')[0],':', unique=False)  # Return values are assumed to all be unique.  Unique to false ensure order is preserved
        self.rank = len(self.LHS)
        
        # Make a lambda for this equation
        Srhs = fortran3d( self.eqstr.split('=')[1] , sMap)
        self.RHS = eval( 'lambda self: ' + Srhs )

        
        # Check to see if this is conserved PDE
        if ( 'ddt(' in eqstr ):
            self.kind = 'PDE'
        
            
        
                     
                
            
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
            

    def emptyScalar(self):
        return numpy.empty( self.chunk_3d_size, dtype=numpy.float64, order='F')*0.0

    
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
    
    
class pyrandaSim:

    def __init__(self,name,meshOptions):

        self.name = name

        self.meshOptions = meshOptions

        if meshOptions.has_key('type'):
            self.meshType = meshOptions['type']
        else:            
            raise ValueError('No suitable mesh type specified.')

        if self.meshType == 'cartesian':

            try:
                nx = meshOptions['nn'][0]
                ny = meshOptions['nn'][1]
                nz = meshOptions['nn'][2]

                dx = (meshOptions['xn'][0]-meshOptions['x1'][0])/max(nx-1,1)
                dy = (meshOptions['xn'][1]-meshOptions['x1'][1])/max(ny-1,1)
                dz = (meshOptions['xn'][2]-meshOptions['x1'][2])/max(nz-1,1)
                periodic = meshOptions['periodic']

                self.PyMPI = pyrandaMPI(nx,ny,nz,dx,dy,dz,periodic)

            except:
                raise ValueError("Invalid options given for cartesian mesh")

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.mesh = pyrandaMesh(self.meshOptions,self.PyMPI)
        
        self.equations = []
        self.conserved = []
        self.variables = {}
        self.nPDE = 0
        self.nALG = 0

        # Package info
        self.packages = {}

        # Compute sMap
        self.get_sMap()

        # Time
        self.time = 0.0

    def iprint(self,sprnt):
        if self.PyMPI.master:
            print(sprnt)
        
    def addPackage(self,package):
        
        self.packages[package.name] = package
        package.get_sMap()
        self.sMap.update( package.sMap )
            
    def allocate(self):
        """
        Loop over vars and allocate the data
        """
        for v in self.variables:
            self.variables[v].allocate(self.PyMPI)
        

    def addVar(self,name,kind=None,rank='scalar'):
        var = pyrandaVar(name,kind,rank)
        self.variables[name] = var

    def addEqu(self,equation):

        peq = pyrandaEq(equation,self.sMap) 
        self.equations.append( peq )
        if peq.kind == 'PDE':
            self.nPDE += 1
            self.conserved.append( peq.LHS[0] )
            if len( peq.LHS ) > 1:
                print 'Warning... only single return values for PDEs allowed'
                exit()
        elif peq.kind == 'ALG':
            self.nALG += 1
        else:
            raise ValueError('No suitable equation type specified for %s' % peq.eqstr)

    def EOM(self,eom):
        """
        Higher level wrapper to make equations of motion from a single string
        Will add eqautions and variables as needed.
        """

        # Split up the equation lines
        eom_lines = filter(None,eom.split('\n'))
        eom_lines = [el for el in eom_lines if el.strip()[0] != '#']  # Comments work
        var_names = []

        #### Get unique set of variables ####
        for eq in eom_lines:
            evars = findVar( eq, ':')
            var_names += evars

        var_names = list(set(var_names))

        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classisfy variables
        
        if self.PyMPI.master:
            for variable in self.variables:
                print 'Adding variables: ' , variable

        #### Add equations of motion ####
        if self.PyMPI.master:
            print 'Adding equations of motion: '

        eqN = 1
        for eq in eom_lines:
            self.addEqu( eq )
            if self.PyMPI.master:
                print '(%s)   %s' % (eqN,eq)
            eqN += 1
            

        # Set up the variables in memory
        self.allocate()
            

        
    def emptyScalar(self):
        return self.PyMPI.emptyScalar()
        
    def updateFlux(self):
        #
        ncons = self.nPDE
        shape = self.mesh.shape[:]
        shape.append( ncons )        
        # Only ddt() terms are extracted
        flux = {} #numpy.asfortranarray( numpy.zeros( shape ) )
        #ieq = 0
        for eqo in self.equations:
            eq = eqo.eqstr
            if ( eqo.kind == 'PDE' ):
                lhs = eqo.LHS[0]
                Srhs = eq.split('=')[1]  # This is a string to evaluate
                if USE_LAMBDA:
                    flux[lhs] = eqo.RHS(self)
                else:
                    rhs = eval( fortran3d(Srhs,self.sMap) )  # This is where the work is done
                    flux[lhs] = rhs

        return flux

    def updateVars(self):
        #
        # Update the equations
        #
        for eq in self.equations:
            if ( eq.kind == 'ALG'):
                if True:
                    rhs = eq.RHS(self)
                    if eq.rank == 1:
                        self.variables[eq.LHS[0]].data = rhs
                    else:
                        for ii in range(len(rhs)):
                            #import pdb
                            #pdb.set_trace()
                            self.variables[eq.LHS[ii]].data = rhs[ii]
                else:
                    exec( fortran3d(eq.eqstr,self.sMap) )
            
        
                
        
    def ddx(self,val):
        if self.nx <= 1:
            return 0.0
        dfdx = self.emptyScalar()
        try:
            self.PyMPI.der.ddx( val, dfdx )
        except:
            import pdb
            pdb.set_trace()

        return dfdx

    def ddy(self,val):
        if self.ny <= 1:
            return 0.0
        dfdy = self.emptyScalar()
        self.PyMPI.der.ddy( val, dfdy )
        return dfdy

    def ddz(self,val):
        if self.nz <= 1:
            return 0.0
        dfdz = self.emptyScalar()
        self.PyMPI.der.ddz( val, dfdz )
        return dfdz
    
    def grad(self,val):
        return [self.ddx(val),self.ddy(val),self.ddz(val)]

    def laplacian(self,val):
        return self.PyMPI.der.laplacian( val )

    def filterx(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_x(val, f_tilde)
        return f_tilde

    def filtery(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_y(val, f_tilde)
        return f_tilde

    def filterz(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_z(val, f_tilde)
        return f_tilde 

    def filter(self,val):
        if self.nx > 1:
            f1 = self.filterx(val)
        else:
            f1 = val
        if self.ny > 1:
            f2 = self.filtery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.filterz(f2)
        else:
            f1 = f2
        return f1

    def gfilterx(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_x(val, f_tilde)
        return f_tilde

    def gfiltery(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_y(val, f_tilde)
        return f_tilde

    def gfilterz(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_z(val, f_tilde)
        return f_tilde 

    def gfilter(self,val):
        if self.nx > 1:
            f1 = self.gfilterx(val)
        else:
            f1 = val
        if self.ny > 1:            
            f2 = self.gfiltery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.gfilterz(f2)
        else:
            f1 = f2
        return f1

    
    def rk4(self,time,dt):

        Ark = [0.0]*5
        Ark[0] = 0.0;
        Ark[1] = -6234157559845./12983515589748.;
        Ark[2] = -6194124222391./4410992767914.;
        Ark[3] = -31623096876824./15682348800105.;
        Ark[4] = -12251185447671./11596622555746.;

        Brk = [0.0]*5
        Brk[0] = 494393426753./4806282396855.;
        Brk[1] = 4047970641027./5463924506627.;
        Brk[2] = 9795748752853./13190207949281.;
        Brk[3] = 4009051133189./8539092990294.;
        Brk[4] = 1348533437543./7166442652324.;

        eta = [0.0]*5
        eta[0] = 494393426753./4806282396855.;
        eta[1] = 4702696611523./9636871101405.;
        eta[2] = 3614488396635./5249666457482.;
        eta[3] = 9766892798963./10823461281321.;
        eta[4] = 1.0;

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        tmp1 = {}
        tmp2 = {}
        PHI  = {}
        for U in self.conserved: 
            tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        # Get primative flow variables
        self.updateVars()
        time_i = time
        for ii in range(5):
            #    ii
            FLUX = self.updateFlux()
            for U in self.conserved:
                tmp1[U] =  Ark[ii]*PHI[U]
                PHI[U]  =  dt*FLUX[U] + tmp1[U]
                tmp2[U] =  Brk[ii]*PHI[U]
                self.variables[U].data =  self.variables[U].data + tmp2[U]
            time = time_i + eta[ii]*dt
            self.time = time
            self.updateVars()
        return time

    def get_sMap(self):
        sMap = {}
        sMap['ddx(' ] = 'self.ddx('
        sMap['ddy(' ] = 'self.ddy('
        sMap['ddz(' ] = 'self.ddz('
        sMap['fbar('] = 'self.filter('
        sMap['gbar('] = 'self.gfilter('
        sMap['grad('] = 'self.grad('
        sMap['simtime'] = 'self.time'
        sMap['dot(' ] = 'numpy.dot('
        sMap['abs(' ] = 'numpy.abs('
        sMap['sqrt(' ] = 'numpy.sqrt('
        sMap['sin('] = 'numpy.sin('
        sMap['tanh('] = 'numpy.tanh(' 
        sMap['lap(' ] = 'self.laplacian('
        
        sMap[':x:']   = 'self.mesh.coords[0]'
        sMap[':y:']   = 'self.mesh.coords[1]'
        sMap[':z:']   = 'self.mesh.coords[2]'
        self.sMap = sMap
        
    
            
def fortran3d(form,sMap):
    """
    Formula translator: Take in a string and return 
    the mapped string for evaluation
    """
    keyS = ':'
   
    # Fill in form with variable data arrays
    varList = findVar(form,keyS)
    for ivar in varList:
        sMap[':%s:'%ivar] = 'self.variables["%s"].data'%ivar
    #
    for mvar in sMap:
        form = form.replace(mvar,sMap[mvar])
    #
    return form

    
    
def findVar(string,keyS,unique=True):
    """
    Return a list of the variable names to replace
    """
    # check that there are even numbers
    cnt = string.count(keyS)
    #
    if (cnt % 2) != 0:
        print "Error: unclosed variable definition marker ::"
        raise ValueError(string)
    #
    varStr = []
    s = ''
    inc = False
    for ss in string:
        if ss == keyS:
            inc = not inc
            if not inc:
                varStr.append(s)
                s = ''
            else:
                continue
        if inc:
            s += ss
    #
    # list(set makes LHS our of order
    if unique:
        varStr = list(set(varStr))
    
    return varStr

