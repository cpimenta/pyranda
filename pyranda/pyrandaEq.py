from pyrandaUtils import *

        
class pyrandaEq:

    def __init__(self,eqstr,sMap,pympi):
        """
        Read in eqstr and extract the LHS variable 
        and the kind of equation (PDE or ALG)
        """
        self.eqstr = eqstr
        self.kind = 'ALG'
        self.active = True
        self.rank = 1
        self.sRHS = ''
        if '=' in eqstr:
            self.LHS = findVar(eqstr.split('=')[0],'scalar', unique=False)  # Return values are assumed to all be unique.  Unique to false ensure order is preserved
            self.rank = len(self.LHS)
        else:
            self.LHS = None   # No return functions
            
                       
        # Make a lambda for this equation
        if self.LHS:
            Srhs = fortran3d( self.eqstr.split('=')[1] , sMap)
        else:
            Srhs = fortran3d( self.eqstr , sMap)
        self.sRHS = Srhs
        self.RHS = eval( 'lambda self: ' + Srhs )

        
        # Check to see if this is conserved PDE
        if ( 'ddt(' in eqstr ):
            self.kind = 'PDE'
