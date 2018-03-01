import numpy


class pyrandaVar:

    def __init__(self,name,kind,rank):

        self.name = name
        self.kind = kind
        self.rank = rank

        self.rhs = None   # RHS for algebraic closure
        self.dt  = None   # dt for temporal integration
        self.allocated = False
        
    def allocate(self,pympi):


        # Only allow single allocation
        if not self.allocated:
            self.pydata = pyrandaData(self.rank,pympi)
        
            # Inherit the mesh size
            #if self.rank == 'scalar':
            #    self.data = pympi.emptyScalar()
            #elif self.rank == 'vector':
            #    self.data = pympi.emptyVector()
            #elif self.rank == 'tensor':
            #    self.data = pympi.emptyTensor()
            #else:
            #    raise ValueError('Error: rank: %s not recognized' % self.rank)

        self.allocated = True


class pyrandaData:

    def __init__(self,rank,pympi):

        self.pympi = pympi

        # Convert string rank to integer
        if type(rank) == type('rank'):
            if rank == 'scalar':
                self.rank = 1
            elif rank == 'vector':
                self.rank = 2
            elif rank == 'tensor':
                self.rank = 3
        else:
            self.rank = rank
            
        # Inherit the mesh size
        if self.rank == 1:
            self.data = [pympi.emptyScalar()]
        elif self.rank == 2:
            self.data = pympi.emptyVector()
        elif self.rank == 3:
            self.data = pympi.emptyTensor()
        else:
            raise ValueError('Error: rank: %s not recognized' % self.rank)



    def get(self):
        return self.data

    def set(self,sdata):
        self.data = sdata

    # Simple math function overloads


    def __mul__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = self.data[0] * other.data[rr]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = other.data[0] * self.data[rr]
                return tmp
            else:                
                print "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] * other
            return tmp

                
        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] * other
            return tmp

        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'").format(self.__class__, type(other))
            


    __rmul__ = __mul__
    
    
    def __div__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = self.data[0] / other.data[rr]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = self.data[rr] / other.data[0]
                return tmp
            else:
                string =  "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)      
        

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] / other
            return tmp
        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] / other
            return tmp
        
    def __rdiv__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = other.data[rr] / self.data[0]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = other.data[0] / self.data[rr]
                return tmp
            else:
                string =  "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)      
        

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = other / self.data[rr]
            return tmp

        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = other / self.data[rr]
            return tmp


    def __add__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = self.data[0] + other.data[rr]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = self.data[rr] + other.data[0]
                return tmp
            else:
                string =  "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)      
        

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] + other
            return tmp

        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] + other
            return tmp

    __radd__ = __add__
        


    def __sub__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = self.data[0] - other.data[rr]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = self.data[rr] - other.data[0]
                return tmp
            else:
                string =  "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)      
        

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] - other
            return tmp

        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = self.data[rr] - other
            return tmp


    def __rsub__(self,other):

        if isinstance(other,self.__class__):
            if self.rank == 1:
                tmp = pyrandaData(other.rank,other.pympi)
                for rr in range(other.rank):
                    tmp.data[rr] = other.data[rr] - self.data[0]
                return tmp
            elif other.rank == 1:
                tmp = pyrandaData(self.rank,self.pympi)
                for rr in range(self.rank):
                    tmp.data[rr] = other.data[0] - self.data[rr]
                return tmp
            else:
                string =  "Error: Scalar-Vector,Scalar-Tensor,only allowed"
                raise ValueError(string)      
        

        elif isinstance(other,int):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = other - self.data[rr]
            return tmp

        elif isinstance(other,float):
            tmp = pyrandaData(self.rank,self.pympi)
            for rr in range(self.rank):
                tmp.data[rr] = other - self.data[rr]
            return tmp


    def __neg__(self):


        tmp = pyrandaData(self.rank,self.pympi)

        for rr in range(self.rank):
            tmp.data[rr] = -self.data[0] 
        return tmp


    def __pow__(self,a):

        tmp = pyrandaData(self.rank,self.pympi)

        for rr in range(self.rank):
            tmp.data[rr] = self.data[0] ** a
        return tmp


    def __lt__(self,other):

        tmp = pyrandaData(self.rank,self.pympi)
        
        if type(other) == type(1.0):
            oth = other
        else:
            oth = other.data[0]
        
        for rr in range(self.rank):
            tmp.data[rr] = self.data[0] < oth
        return tmp

    def __gt__(self,other):

        tmp = pyrandaData(self.rank,self.pympi)

        for rr in range(self.rank):
            tmp.data[rr] = self.data[0] > other.data[0]
        return tmp


    def __getitem__(self,args):
        import pdb
        pdb.set_trace()
    

    
def pyData_exp(pdata):

    tmp = pyrandaData(pdata.rank,pdata.pympi)
    
    for rr in range(pdata.rank):
        tmp.data[rr] = numpy.exp(pdata.data[0])
    return tmp


def pyData_sqrt(pdata):

    tmp = pyrandaData(pdata.rank,pdata.pympi)
    
    for rr in range(pdata.rank):
        tmp.data[rr] = numpy.sqrt(pdata.data[0])
    return tmp

def pyData_abs(pdata):

    tmp = pyrandaData(pdata.rank,pdata.pympi)
    
    for rr in range(pdata.rank):
        tmp.data[rr] = numpy.abs(pdata.data[0])
    return tmp


def pyData_where(cond,val_if,val_el):

    tmp = pyrandaData(cond.rank,cond.pympi)

    if cond.rank != 1:
        print "Where statement only supports rank 1"


    if type(val_if) == type(1.0):
        ifval = val_if
    else:
        ifval = val_if.data[0]

    if type(val_el) == type(1.0):
        elval = val_el
    else:
        elval = val_el.data[0]

        
    tmp.data[0] = numpy.where( cond.data[0], ifval, elval )

    return tmp
