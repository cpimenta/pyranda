
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
        elif self.rank == 'vector':
            self.data = pympi.emptyVector()
        else:
            raise ValueError('Error: rank: %s not recognized' % self.rank)
