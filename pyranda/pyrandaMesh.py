import numpy 


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

