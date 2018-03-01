from mpi4py import MPI
import numpy 
import re
import sys
import time
sys.path.append('/Users/olson45/Research/FloATPy')

from floatpy.parallel import t3dmod
from floatpy.derivatives.compact import CompactDerivative
from floatpy.filters.filter import Filter

# Import the pyranda library
from pyrandaMPI  import pyrandaMPI
from pyrandaVar   import *
from pyrandaEq    import pyrandaEq
from pyrandaMesh  import pyrandaMesh
from pyrandaUtils import *

            
                                              
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
        self.initial_conditions = []
        self.nPDE = 0
        self.nALG = 0
        self.nIC  = 0

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
        if name not in self.variables:
            var = pyrandaVar(name,kind,rank)
            self.variables[name] = var

    def addEqu(self,equation,IC=False):

        peq = pyrandaEq(equation,self.sMap,self.PyMPI)

        # Check to see if this is conserved PDE
        if ( 'ddt(' in equation ):
            peq.kind = 'PDE'

        # Check for IC flag
        if IC:
            peq.kind = 'IC'
            
        self.equations.append( peq )
        if peq.kind == 'PDE':
            self.nPDE += 1
            self.conserved.append( peq.LHS[0] )
            if len( peq.LHS ) > 1:
                print 'Warning... only single return values for PDEs allowed'
                exit()
        elif peq.kind == 'ALG':
            self.nALG += 1
        elif peq.kind == 'IC':
            self.nIC  += 1
        else:
            raise ValueError('No suitable equation type specified for %s' % peq.eqstr)

    def EOM(self,eom):
        """
        Higher level wrapper to make equations of motion from a single string
        Will add eqautions and variables as needed.
        """
        vKey = ':'
        
        # Split up the equation lines
        eom_lines = filter(None,eom.split('\n'))
        eom_lines = [el.replace(' ','') for el in eom_lines ]

        eom_lines = filter(None,eom_lines)
        eom_lines = [el for el in eom_lines if el.strip()[0] != '#']  # Comments work
        var_names = []

        #### Get unique set of variables ####
        for eq in eom_lines:
            evars = findVar( eq, vKey)
            var_names += evars

        var_names = list(set(var_names))

        # Of the var_names, detect vectors and tensors (default is scalar)
        var_names_vector = []
        var_names_tensor = []
        for eq in eom_lines:
            for iiv in var_names:
                tmp = eq.replace('%s%s%s'%(vKey,iiv,vKey),'').strip()
                if tmp == 'vec()':
                    var_names.remove(iiv)
                    var_names_vector.append( iiv )
                if tmp == 'ten()':
                    var_names.remove(iiv)
                    var_names_tensor.append( iiv )

        for svar in var_names:
            self.addVar(svar,kind='conserved')   # Todo: classify variables

        for vvar in var_names_vector:
            self.addVar(vvar,kind='conserved',rank='vector')   

        for tvar in var_names_vector:
            self.addVar(tvar,kind='conserved',rank='tensor')   
            
        if self.PyMPI.master:
            for variable in self.variables:
                print 'Adding variables: ' , variable

        #### Add equations of motion ####
        if self.PyMPI.master:
            print 'Adding equations of motion: '

        eqN = 1
        for eq in eom_lines:
            # Filter out the vec/tensor declarations
            tmp = ''
            try:
                tmp = eq.split('(')[0].strip()
            except:
                pass
            if ('vec' in tmp) or ('ten' in tmp):
                continue
            
            self.addEqu( eq )
            if self.PyMPI.master:
                print '(%s)   %s' % (eqN,eq)
            eqN += 1
            

        # Set up the variables in memory
        self.allocate()
            
    def setIC(self,ics):
        """
        Evaluate the initial conditions and then update variables
        """
        # Split up the equation lines
        ic_lines = filter(None,ics.split('\n'))
        ic_lines = [el.replace(' ','') for el in ic_lines ]  # Comments work
        ic_lines = filter(None,ic_lines)
        ic_lines = [el for el in ic_lines if el.strip()[0] != '#']  # Comments work

        var_names = []

        #### Get unique set of variables ####
        for eq in ic_lines:
            evars = findVar( eq, ':')
            var_names += evars

        var_names = list(set(var_names))

        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classify variables
        
        for ic in ic_lines:
            ic_mod = ic #+ '+self.emptyScalar()'
            self.addEqu(ic,IC=True)
            #exec(fortran3d(ic_mod,self.sMap))                           

        self.allocate()
        self.updateVars(kind="IC")    # Enforce the IC
        self.updateVars()             # Update algebraic closures
        
        
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
                try:
                    lhs = eqo.LHS[0]
                    Srhs = eq.split('=')[1]  # This is a string to evaluate
                    flux[lhs] = eqo.RHS(self)
                except:
                    import pdb
                    pdb.set_trace()

        return flux

    def updateVars(self,kind='ALG'):
        #
        # Update the algebraic equations
        #
        for eq in self.equations:

            if not eq.active:
                continue

            if ( eq.kind == kind):
                try:
                    rhs = eq.RHS(self)

                    # If no LHS , then on to the next eq
                    if not eq.LHS:
                        continue

                    # Only rank one updates for now...
                    if eq.rank == 1:
                        self.variables[eq.LHS[0]].pydata = rhs
                        
                except:
                    print("Error updating equation: %s" % eq.eqstr)
                    import pdb
                    pdb.set_trace()
            
                    
    def ddx(self,val):

        pyd = False   # If val is pyrandaVar object

        if self.nx <= 1:
            return 0.0

        dfdx = self.emptyScalar()
        
        try:
            if val.__module__ == 'pyranda.pyrandaVar':
                self.PyMPI.der.ddx( val.data[0], dfdx )
                pyd = True
            else:
                self.PyMPI.der.ddx( val, dfdx )
            
        except:
            import pdb
            pdb.set_trace()


        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = dfdx
            return tmp
        else:
            return dfdx

    def ddy(self,val):
        pyd = False
        if self.ny <= 1:
            return 0.0
        dfdy = self.emptyScalar()
        if val.__module__ == 'pyranda.pyrandaVar':
            self.PyMPI.der.ddy( val.data[0], dfdy )
            pyd = True
        else:
            self.PyMPI.der.ddy( val, dfdy )

        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = dfdy
            return tmp
        else:
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

        pyd = False
        if val.__class__.__name__ == "pyrandaData":
            pyd = True
            tmp = pyrandaData(val.rank,val.pympi)
            lap = self.PyMPI.der.laplacian( val.data[0] )
            tmp.data[0] = lap
            return tmp
        else:        
            return self.PyMPI.der.laplacian( val )

    def filterx(self,val):
        pyd = False
        f_tilde = self.emptyScalar()
        if val.__class__.__name__ == 'pyrandaData':
            self.PyMPI.fil.filter_x(val.data[0], f_tilde)
            pyd = True
        else:
            self.PyMPI.fil.filter_x(val, f_tilde)
        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = f_tilde
            return tmp
        else:
            return f_tilde

    def filtery(self,val):
        pyd = False
        f_tilde = self.emptyScalar()
        if val.__class__.__name__ == 'pyrandaData':
            self.PyMPI.fil.filter_y(val.data[0], f_tilde)
            pyd = True
        else:
            self.PyMPI.fil.filter_y(val, f_tilde)

        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = f_tilde
            return tmp
        else:
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
        pyd = False
        f_tilde = self.emptyScalar()
        if val.__module__ == 'pyranda.pyrandaVar':
            self.PyMPI.gfil.filter_x(val.data[0], f_tilde)
            pyd = True
        else:
            self.PyMPI.gfil.filter_x(val, f_tilde)
        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = f_tilde
            return tmp
        else:
            return f_tilde

    def gfiltery(self,val):
        pyd = False
        f_tilde = self.emptyScalar()
        if val.__module__ == 'pyranda.pyrandaVar':
            self.PyMPI.gfil.filter_y(val.data[0], f_tilde)
            pyd = True
        else:
            self.PyMPI.gfil.filter_y(val, f_tilde)
        if pyd:
            tmp = pyrandaData(val.rank,val.pympi)
            tmp.data[0] = f_tilde
            return tmp
        else:
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

    def var(self,name,index=0):
        """
        Quick access to mesh variables
        name - name of the variable
        index (optional) index into the rank list
        """
        if name in self.variables.keys():            
            return self.variables[name].pydata.data[index]
        else:
            msg = "%s is not a valid variable.\n" % name
            self.iprint( msg )
        
    
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
            #tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            #tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            #PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp1[U] = pyrandaData(1,self.PyMPI )
            tmp2[U] = pyrandaData(1,self.PyMPI)
            PHI[U]  = pyrandaData(1,self.PyMPI)

            
        # Get primative flow variables
        #self.updateVars()
        time_i = time
        for ii in range(5):
            #    ii
            FLUX = self.updateFlux()
            for U in self.conserved:
                tmp1[U] =  Ark[ii]*PHI[U]
                PHI[U]  =  dt*FLUX[U] + tmp1[U]
                tmp2[U] =  Brk[ii]*PHI[U]
                #self.variables[U].data =  self.variables[U].data + tmp2[U]
                self.variables[U].pydata =  self.variables[U].pydata + tmp2[U]
            time = time_i + eta[ii]*dt
            self.time = time
            self.updateVars()
        return time

    def get_sMap(self):
        sMap = {}
        sMap['ddx(' ] = 'self.ddx('
        sMap['ddxPD(' ] = 'self.ddx('
        sMap['ddy(' ] = 'self.ddy('
        sMap['ddz(' ] = 'self.ddz('
        sMap['fbar('] = 'self.filter('
        sMap['gbar('] = 'self.gfilter('
        sMap['grad('] = 'self.grad('
        sMap['simtime'] = 'self.time'
        sMap['lap(' ] = 'self.laplacian('
        sMap['sum(' ] = 'self.PyMPI.sum3D('
        sMap['sign(' ] = 'numpy.sign('
        sMap['dot(' ] = 'numpy.dot('
        sMap['abs(' ] = 'pyData_abs('
        sMap['sqrt(' ] = 'pyData_sqrt('
        sMap['sin('] = 'numpy.sin('
        sMap['tanh('] = 'numpy.tanh('
        sMap['exp('] = 'pyData_exp(' 
        sMap['where('] = 'pyData_where('
        sMap['3d()'] = 'self.emptyScalar()'
        sMap[':pi:'] = 'numpy.pi'
        
        sMap[':x:']   = 'self.mesh.coords[0]'
        sMap[':y:']   = 'self.mesh.coords[1]'
        sMap[':z:']   = 'self.mesh.coords[2]'
        self.sMap = sMap
        


    def zbar(self,var):

        pyd = self.variables[var].pydata
        if pyd.rank == 1:
            return self.PyMPI.zbar( pyd.data[0] )
        else:
            val = []
            for rr in prd.rank:
                val.append( self.pyMPI.zbar( pyd.data[rr] ) )
            return val
                            
