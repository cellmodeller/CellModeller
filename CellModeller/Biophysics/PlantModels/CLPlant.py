import pyopencl as cl
import pyopencl.array as cl_array
import numpy
import numpy.linalg as la
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel
import math
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.raw.GL.VERSION.GL_1_5 import glBufferData as rawGlBufferData

class CLPlant:
    """ Implementation of CM3 model using OpenCL and matrix free FEM
    """

    maxNodes = 16384*4 
    maxCells =16384
    maxWalls = 16384*4
    maxNbs = numpy.uint32(16)
    maxCellWalls = 64 
    maxCellNodes = 64 

    nextNodeIdx = 0
    nextWallIdx = 0
    nextCellIdx = 0

    Karea = 0.1
    Kperim = 1.0

    # Units for params:
    # S = wall cross-section (um^2)
    # I = moment of inertia (um^4)
    # P = turgor pressure (N/um^2)
    # mu = viscosity (N s/um^2)
    # ---

    def __init__(self, tikparam, cgstol, S=5.0, I=2e2, P=0.5e-6, mu=12.0, useVBO=False):
        # FEM params
        self.S=numpy.float32(S)
        self.I=numpy.float32(I)
        self.P=numpy.float32(P*40.0) # multiply by estimate of cell height
        self.mu=numpy.float32(mu)

        self.useVBO = useVBO
        
        self.initOpenCLContext()
        # arrays
        self.nodes = numpy.zeros((self.maxNodes,)).astype(cl_array.vec.float4)
        self.disp = numpy.zeros((self.maxNodes,)).astype(cl_array.vec.float4)
        self.cell_nodes =  numpy.zeros((self.maxCells,self.maxCellNodes+1)).astype(numpy.int32)
        self.cell_nodes[:,0] = -1 # no cells
        self.cell_walls=  numpy.zeros((self.maxCells,self.maxCellWalls+1)).astype(numpy.int32)
        self.cell_walls[:,0] = -1 # no cells
        self.cell_turgor=  numpy.zeros((self.maxCells,)).astype(numpy.float32)
        self.wall_nodes = numpy.zeros((self.maxWalls,2)).astype(numpy.uint32)
        self.wall_cells = numpy.zeros((self.maxWalls,2)).astype(numpy.int32)
        self.wall_cells.fill(-1) # no walls initially
        self.wall_stiff = numpy.zeros((self.maxWalls,)).astype(numpy.float32)
        self.node_nbs = numpy.zeros((self.maxNodes,self.maxNbs+1)).astype(numpy.int32)
        self.node_walls = numpy.zeros((self.maxNodes,self.maxNbs+1)).astype(numpy.int32)
        self.node_walls[:,0]=-1 # no walls initially

        self.talpha = numpy.float32(tikparam)
        self.cgsTol = cgstol
        self.createKernels()
        #self.createGapCells()
        #self.create1Cell(10.0, 4)
        self.create2Cells()
        #self.createBeam(1,1,1)
        self.initDeviceData()

    def initOpenCLContext(self):
        # NB. Seems there needs to be a current GL context before setting up gl 
        # interop here
        from pyopencl.tools import get_gl_sharing_context_properties
        # Some OSs prefer clCreateContextFromType, some prefer
        # clCreateContext. Try both.
        platform = cl.get_platforms()[0]
        try:
            self.ctx = cl.Context(properties=[
            (cl.context_properties.PLATFORM, platform)]
            + get_gl_sharing_context_properties())
        except:
            self.ctx = cl.Context(properties=[
            (cl.context_properties.PLATFORM, platform)]
            + get_gl_sharing_context_properties(),
            devices = [platform.get_devices()[0]])
        #self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

        mf = cl.mem_flags

    def createKernels(self):
        # The source for kernels that compute Ax and b, in solving Ax=b
        kernelSource = open('CLPlant.cl', 'r').read()
        self.FEMKernels = cl.Program(self.ctx, kernelSource).build()
        # Some kernels that seem like they should be built into pyopencl...
        self.vclear = ElementwiseKernel(self.ctx, "float4 *v", "v[i]=0.0", "vecclear")
        self.vadd = ElementwiseKernel(self.ctx, "float4 *res, const float4 *in1, const float4 *in2", \
                "res[i] = in1[i] + in2[i]", "vecadd")
        self.vsub = ElementwiseKernel(self.ctx, "float4 *res, const float4 *in1, const float4 *in2", \
                "res[i] = in1[i] - in2[i]", "vecsub")
        self.vaddkx = ElementwiseKernel(self.ctx, \
                "float4 *res, const float k, const float4 *in1, const float4 *in2", \
                "res[i] = in1[i] + k*in2[i]", "vecaddkx")
        self.vsubkx = ElementwiseKernel(self.ctx, \
                "float4 *res, const float k, const float4 *in1, const float4 *in2", \
                "res[i] = in1[i] - k*in2[i]", "vecsubkx")
        self.vaddfloat4 = ElementwiseKernel(self.ctx, \
                "float4 *res, const float4 k, const float4 *in1", \
                "res[i] = in1[i] + k", "vecaddfloat4")

        # A dot product as sum of float4 dot products - 
        # i.e. like flattening vectors of float4s into big float vectors
        # then computing dot
        self.vdot = ReductionKernel(self.ctx, numpy.float32, neutral="0",
            reduce_expr="a+b", map_expr="dot(x[i],y[i])",
            arguments="__global float4 *x, __global float4 *y")


    def initDeviceData(self):
        sz = self.nodes.shape
        print "sz"
        print sz
        # Create device storage for intermediate variables
        # node displacements
        self.disp_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)

        # node forces
        self.force_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)

        # CGS temporary storage vars
        # res, rhs, p and Ap storage for cgs
        self.res_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)
        self.rhs_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)
        self.p_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)
        self.Ap_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)

        # Energy min temp storage
        self.energy_deriv_dev = cl_array.zeros(self.queue, sz, cl_array.vec.float4)

        # Wall length storage
        self.wall_length_dev = cl_array.zeros(self.queue, (self.maxWalls,), numpy.float32)
        # Cell geometry
        self.cell_perim_dev = cl_array.zeros(self.queue, (self.maxCells,), numpy.float32)
        self.cell_area_dev = cl_array.zeros(self.queue, (self.maxCells,), numpy.float32)
        self.cell_area0_dev = cl_array.zeros(self.queue, (self.maxCells,), numpy.float32)
        self.cell_centre_dev = cl_array.zeros(self.queue, (self.maxCells,), cl_array.vec.float4)
        self.cell_axis_dev = cl_array.zeros(self.queue, (self.maxCells,), cl_array.vec.float4)

        # Copy host cell geometry etc. to device
        self.nodes_dev = cl_array.to_device(self.queue, self.nodes)
        self.node_nbs_dev = cl_array.to_device(self.queue, self.node_nbs)
        self.node_walls_dev = cl_array.to_device(self.queue, self.node_walls)
        self.wall_cells_dev = cl_array.to_device(self.queue, self.wall_cells)
        self.wall_nodes_dev = cl_array.to_device(self.queue, self.wall_nodes)
        self.wall_stiff_dev = cl_array.to_device(self.queue, self.wall_stiff)
        self.cell_walls_dev = cl_array.to_device(self.queue, self.cell_walls)
        self.cell_turgor_dev = cl_array.to_device(self.queue, self.cell_turgor)
        self.cell_nodes_dev = cl_array.to_device(self.queue, self.cell_nodes)

        if self.useVBO:
            self.vbo = glGenBuffers(1)
            glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
            rawGlBufferData(GL_ARRAY_BUFFER, self.maxNodes * 4 * 4, None, GL_STATIC_DRAW)
            glEnableClientState(GL_VERTEX_ARRAY)
            glVertexPointer(3, GL_FLOAT, 4*4, None)
            self.coords_dev = cl.GLBuffer(self.ctx, cl.mem_flags.READ_WRITE, int(self.vbo))

            self.ibo = glGenBuffers(1)
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ibo)
            rawGlBufferData(GL_ELEMENT_ARRAY_BUFFER, self.maxWalls * 2 * 4, None, GL_STATIC_DRAW)
            self.idx_dev = cl.GLBuffer(self.ctx, cl.mem_flags.READ_WRITE, int(self.ibo))

        # Link an OpenGL vertex buffer to nodes
        """self.vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        rawGlBufferData(GL_ARRAY_BUFFER, self.maxNodes * 4 * 4, None, GL_DYNAMIC_DRAW)

        glEnableClientState(GL_VERTEX_ARRAY)
        glVertexPointer(4, GL_FLOAT, 0, None)
        self.nodes_dev = cl.GLBuffer(self.ctx, cl.mem_flags.READ_WRITE, int(self.vbo))
        cl.enqueue_acquire_gl_objects(self.queue, [self.nodes_dev])
        cl.enqueue_copy(self.queue, self.nodes_dev, self.nodes[0:self.nextNodeIdx])
        self.queue.finish()
        cl.enqueue_release_gl_objects(self.queue, [self.nodes_dev])

        # Link an OpenGL index array to wall_nodes
        self.ibo = glGenBuffers(1)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ibo)
        rawGlBufferData(GL_ELEMENT_ARRAY_BUFFER, self.maxWalls * 2 * 4, None, GL_STATIC_DRAW)

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ibo)
        self.wall_nodes_dev = cl.GLBuffer(self.ctx, cl.mem_flags.READ_WRITE, int(self.ibo))
        cl.enqueue_acquire_gl_objects(self.queue, [self.wall_nodes_dev])
        cl.enqueue_copy(self.queue, self.wall_nodes_dev, self.wall_nodes[0:self.nextWallIdx])
        self.queue.finish()
        cl.enqueue_release_gl_objects(self.queue, [self.wall_nodes_dev])
"""

    def createBeam(self, L, E, P):
        self.nodes[0] = (0,0,0,0)
        self.nodes[1] = (L,0,0,0)
        self.nextNodeIdx = 2

        self.cell_nodes[0,0:3] = [0,1,-1]
        self.cell_walls[0,0:2] = [0,-1]

        self.cell_turgor[0] = P
        self.nextCellIdx = 1

        self.wall_nodes[0] = [0,1]
        self.wall_cells[0] = [0,-1]
        self.wall_stiff[0] = E
        self.nextWallIdx = 1

        self.node_nbs[0:2,0:2] = [[1,-1],[0,-1]]
        self.node_walls[0:2,0:2] = [[0,-1],[0,-1]]

    def createGapCells(self):
        # Set up 2 adjacent square cells with 4 walls in host arrays

        # Node positions
        self.nodes[0] = (0,-100,0,0)
        self.nodes[1] = (100,-100,0,0)
        self.nodes[2] = (100,0,0,0)
        self.nodes[3] = (100,100,0,0)
        self.nodes[4] = (0,100,0,0)
        self.nodes[5] = (-100,100,0,0)
        self.nodes[6] = (-100,0,0,0)
        self.nodes[7] = (-100,-100,0,0)
        self.nodes[8] = (0,-25,0,0)
        self.nodes[9] = (25,0,0,0)
        self.nodes[10] = (0,25,0,0)
        self.nodes[11] = (-25,0,0,0)
        self.nextNodeIdx = 12 


        # Indices of cell nodes
        self.cell_nodes[0:4,0:6] =  [[0,8,11,6,7,-1],[0,1,2,9,8,-1],[2,3,4,10,9,-1],[4,5,6,11,10,-1]]
        self.cell_walls[0:4,0:6] =  [[7,8,15,11,6,-1],[0,1,9,12,8,-1],[9,2,3,10,13,-1],[10,4,5,11,14,-1]]

        # Turgur pressure of cells
        self.cell_turgor[0:4] = self.P

        self.nextCellIdx = 4

        # wall node indices
        self.wall_nodes[0:16] = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,0],[0,8],[2,9],[4,10],[6,11],[8,9],[9,10],[10,11],[11,8]]
        # wall cell indices
        self.wall_cells[0:16] = [[1,-1],[1,-1],[2,-1],[2,-1],[3,-1],[3,-1],[0,-1],[0,-1],[0,1],[1,2],[2,3],[3,0],[-1,1],[-1,2],[-1,3],[-1,0]]
        # wall stiffness
        self.wall_stiff[0:16] = self.mu

        self.nextWallIdx = 16

        # node neighbourhoods
        self.node_nbs[0:12,0:4] = [[1,8,7,-1],[0,2,-1,-1],[3,9,1,-1],[2,4,-1,-1],[3,10,5,-1],[4,6,-1,-1],[5,11,7,-1],[0,6,-1,-1],[0,9,11,-1],[2,8,10,-1],[4,9,11,-1],[6,8,10,-1]]
        # node wall indices
        self.node_walls[0:12,0:4] = [[0,8,7,-1],[0,1,-1,-1],[1,2,9,-1],[2,3,-1,-1],[3,4,10,-1],[4,5,-1,-1],[5,11,6,-1],[6,7,-1,-1],[8,12,15,-1],[9,12,13,-1],[10,13,14,-1],[11,14,15,-1]]


    def create1Cell(self, r, n):
        # Set up one cell with n walls

        # Node positions
        for i in range(n):
            ang = 2.0*i*math.pi/n
            self.nodes[i] = (r*math.cos(ang), r*math.sin(ang), 0, 0)
        self.nextNodeIdx = n


        # Indices of cell nodes
        self.cell_nodes[0,0:n+1] =  range(n) + [-1]
        self.cell_walls[0,0:n+1] =  range(n) + [-1]

        # Turgur pressure of cells
        self.cell_turgor[0] = self.P

        self.nextCellIdx = 1

        # wall node indices
        self.wall_nodes[0:n] = [[i,i+1] for i in range(n-1)] + [[n-1,0]]
        # wall cell indices
        for i in range(n):
            self.wall_cells[i] = [0,-1]
        # wall stiffness
        self.wall_stiff[0:n] = self.mu

        self.nextWallIdx = n

        # node neighbourhoods
        self.node_nbs[0:n,0:3] = [[n-1,1,-1]] + [ [i,i+1,-1] for i in range(n-1) ]
        # node wall indices
        self.node_walls[0:n,0:3] = [[n-1,0,-1]] + [ [i,i+1,-1] for i in range(n-1) ]


    def create2Cells(self):
        # Set up 2 adjacent square cells with 4 walls in host arrays

        # Node positions
        self.nodes[0] = (0,0,0,0)
        self.nodes[1] = (0,10,0,0)
        self.nodes[2] = (-10,10,0,0)
        self.nodes[3] = (-10,0,0,0)
        self.nodes[4] = (10,0,0,0)
        self.nodes[5] = (10,10,0,0)
        self.nextNodeIdx = 6


        # Indices of cell nodes
        self.cell_nodes[0:2,0:5] =  [[0,1,2,3,-1],[0,4,5,1,-1]]
        self.cell_walls[0:2,0:5] =  [[0,1,2,3,-1],[0,4,5,6,-1]]

        # Turgur pressure of cells
        self.cell_turgor[0:2] = self.P

        self.nextCellIdx = 2

        # wall node indices
        self.wall_nodes[0:7] = [[0,1],[1,2],[2,3],[3,0],[0,4],[4,5],[5,1]]
        # wall cell indices
        self.wall_cells[0:7] = [[0,1],[0,-1],[0,-1],[0,-1],[1,-1],[1,-1],[1,-1]]
        # wall stiffness
        self.wall_stiff[0:7] = self.mu

        self.nextWallIdx = 7

        # node neighbourhoods
        self.node_nbs[0:6,0:4] = [[1,3,4,-1],[0,2,5,-1],[1,3,-1,-1],[2,0,-1,-1],[0,5,-1,-1],[4,1,-1,-1]]
        # node wall indices
        self.node_walls[0:6,0:4] = [[0,3,4,-1],[0,1,6,-1],[1,2,-1,-1],[2,3,-1,-1],[4,5,-1,-1],[5,6,-1,-1]]


    def intersect(self,widx,c,v):
        # wall intersection copied/modified from CM3.py
        # NB this does not deal with case where line intersects at 
        # a vertex!

        # TJR: Fix to handle vertex intersection:
        # Returns (p, node)
        # p=[x,y,z], node=None --> intersection at point p
        # p=None, node={0,1} --> intersection at node {0,1} of wall

        """
        intersect wall with line through point c along vector v 
        return intersection point, or None if fails
        
        see http://exaflop.org/docs/cgafaq/cga1.html    
        AB=A+r(B-A), r in [0,1]
        CD=C+s(D-C), s in [0,1]
        If AB & CD intersect, then

        A+r(B-A)=C+s(D-C), or

        Ax+r(Bx-Ax)=Cx+s(Dx-Cx)
        Ay+r(By-Ay)=Cy+s(Dy-Cy)  for some r,s in [0,1]
        Solving the above for r and s yields

                (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy)
        r = -----------------------------  (eqn 1)
                (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
                (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------  (eqn 2)
                (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
        Let P be the position vector of the intersection point, then

        P=A+r(B-A) or

        Px=Ax+r(Bx-Ax)
        Py=Ay+r(By-Ay)
        """

        na=self.nodes[self.wall_nodes[widx,0]] 
        a = numpy.array(list(na))[0:2]
        nb=self.nodes[self.wall_nodes[widx,1]] 
        b = numpy.array(list(nb))[0:2]
        da = numpy.array(list(self.disp[self.wall_nodes[widx,0]]))
        db = numpy.array(list(self.disp[self.wall_nodes[widx,1]]))
	L = numpy.linalg.norm(a-b)
        d=c+v

        den=(b[0]-a[0])*(d[1]-c[1])-(b[1]-a[1])*(d[0]-c[0]) #all 2D positions are numpy array shape(2,1), so second [0] needed

        if abs(den)<1e-16: # wall and line are parallel
            print "den too small, %f"%den
            return (None,None,None)

        r=((a[1]-c[1])*(d[0]-c[0])-(a[0]-c[0])*(d[1]-c[1]))/den

        tol=1.0 # the smallest tolerance causes small wall segment which create numerical singularity in the physics
        #print "division r = %e"%r
        #print "division r-1 = %e"%(r-1.0)
        if r<-1e-6 or r>1.0+1e-6:
            return (None,None,None) # No intersection
        elif L*r>tol and L*r<(L-tol):
            # Intersection of wall
            p=a+(b-a)*r
	    d = db*r + (1-r)*da
            return (p, d, None)
        elif L*r<=tol:
            # Intersection at vertex
            return (None, None, 0)
        elif L*r>=(L-tol):
            return (None, None, 1)
        else:
            return (None,None,None)

    def addNodeToCell(self, cidx, nidx):
        cn=0
        while(self.cell_nodes[cidx,cn]>=0):
            cn+=1
        self.cell_nodes[cidx,cn]=nidx
        self.cell_nodes[cidx,cn+1] = -1

    def addWallToCell(self, cidx, widx):
        cw=0
        while(self.cell_walls[cidx,cw]>=0):
            cw+=1
        self.cell_walls[cidx,cw]=widx
        self.cell_walls[cidx,cw+1] = -1

    def whichSideWall(self,widx, pos):
        """
        Determines which side of a wall a point lies by returning 
        the cross product of (p-node0) x (node1 - node0)
        ie +ve one side, -ve the other
        """
        nn0=self.nodes[self.wall_nodes[widx,0]] 
        n0 = numpy.array([nn0[i] for i in range(2)])
        nn1=self.nodes[self.wall_nodes[widx,1]] 
        n1 = numpy.array([nn1[i] for i in range(2)])
        p = numpy.array([pos[i] for i in range(2)])
        #print "n0 = " + str(n0)
        #print "n1 = " + str(n1)
        #print "p = " + str(p)

        #x=(p[1]-n0[1])*(n1[0]-n0[0])-(p[0]-n0[0])*(n1[1]-n0[1])
        #return x[0]
        return numpy.cross((p-n0).transpose(),(n1-n0).transpose())

    def divide_cell(self, cidx, d1id, d2id, centre, axis):
        """
        Divide cell along axis passing through centre
        """
        # find wall intersections
        intersections=[]
        inodes = []
        for w in self.cell_walls[cidx,:]:
            if w==-1:
                break
            v = self.intersect(w,centre,axis)
            #print "v = "+str(v)
            (p, d, nidx)=v
            if p!=None and d!=None:
                intersections.append((w,p,d,nidx))
            if nidx!=None:
                wn = self.wall_nodes[w,nidx]
                # Intersection at a node can occur with 2 walls
                # -- ignore duplicate nodes
                if wn not in inodes:
                    inodes.append(wn)
                    intersections.append((w,p,d,nidx))

        #print "intersections:"
        #print intersections
        if len(intersections)!=2:
            print "intersections:"
            print intersections
            print "axis"
            print axis
            print "centre"
            print centre
            print "nodes"
            print self.nodes[0:self.nextNodeIdx]
            print "walls"
            print self.wall_nodes[0:self.nextWallIdx]
            raise Exception,"Cell division axis does not divide 2 walls"

        # Split intersected walls if needed
        divWallEndNodes=[]
        for w,p,d,nidx in intersections:
            if nidx!=None:
                # division wall hits existing vertex
                divWallEndNodes.append(self.wall_nodes[w,nidx])
            else:
                # Divides wall, split it in 2
                if self.wall_cells[w,0]!=cidx:
                    neighbourCell = self.wall_cells[w,0]
                else: 
                    neighbourCell = self.wall_cells[w,1]
                # add new node and refs
                newNode = self.nextNodeIdx
                self.nextNodeIdx += 1
                self.nodes[newNode] = tuple(p)+(0,0,) # tuple maps to OpenCL float4
		self.disp[newNode] = tuple(d)

                # Add new node to this cell and neighbour
                self.addNodeToCell(cidx, newNode)
                if neighbourCell>=0:
                    self.addNodeToCell(neighbourCell, newNode)
                # Split wall
                oldNode = self.wall_nodes[w,1] 
                self.wall_nodes[w,1] = newNode
                # New wall
                newWall = self.nextWallIdx
                self.nextWallIdx += 1
                self.wall_nodes[newWall] = [newNode, oldNode]
                self.wall_cells[newWall] = self.wall_cells[w]
                #print "wall_cells[newWall] = "
                #print self.wall_cells[newWall]
                self.wall_stiff[newWall] = self.mu
                # Add walls to new node
                self.node_walls[newNode,0:3] = [w, newWall,-1]

                # Replace w with newWall in oldNodes neighbours
                for i in range(self.maxNbs):
                    if self.node_walls[oldNode,i] == w:
                        self.node_walls[oldNode,i] = newWall

                # add new wall
                if neighbourCell>=0:
                    self.addWallToCell(neighbourCell, newWall)
                self.addWallToCell(cidx, newWall)

                divWallEndNodes.append(newNode)

        # New cell
        newCell = self.nextCellIdx
        self.nextCellIdx += 1

        # Create string of wall segments
        daughterwalls0=[]
        daughterwalls1=[]
        dwn1 = divWallEndNodes[0]
        dwn2 = divWallEndNodes[1]
        #print "dwn1 = %i"%dwn1
        #print "nodes[dwn1] = "+str(self.nodes[dwn1])
        pdwn1 = numpy.array([self.nodes[dwn1][i] for i in range(2)])
        pdwn2 = numpy.array([self.nodes[dwn2][i] for i in range(2)])
        nSegs = 1
        for wseg in range(nSegs-1):
            # add new node and refs
            newNode = self.nextNodeIdx
            self.nextNodeIdx += 1
            p = pdwn1 + (pdwn2-pdwn1)*(wseg+1)/nSegs
            self.nodes[newNode] = tuple(p)+(0,) # tuple maps to OpenCL float4
            # Add new node to this cell
            self.addNodeToCell(cidx, newNode)
            wallSeg = self.nextWallIdx
            self.nextWallIdx += 1
            self.wall_nodes[wallSeg] = [dwn1, newNode]
            self.wall_cells[wallSeg] = [cidx, newCell]
            self.wall_stiff[wallSeg] = self.mu
            # Add wall to nodes
            for nw in range(self.maxNbs):
                if self.node_walls[dwn1,nw]==-1:
                    self.node_walls[dwn1,nw] = wallSeg
                    self.node_walls[dwn1,nw+1] = -1
                    break
            for nw in range(self.maxNbs):
                if self.node_walls[newNode,nw]==-1:
                    self.node_walls[newNode,nw] = wallSeg
                    self.node_walls[newNode,nw+1] = -1
                    break
            # Start list of daughter walls with dividing segments
            daughterwalls0.append(wallSeg)
            daughterwalls1.append(wallSeg)
            dwn1 = newNode
        wallSeg = self.nextWallIdx
        self.nextWallIdx += 1
        self.wall_nodes[wallSeg] = [dwn1, dwn2]
        self.wall_cells[wallSeg] = [cidx, newCell] #[newCell,cidx]
        self.wall_stiff[wallSeg] = self.mu
        # Add wall to nodes
        for nw in range(self.maxNbs):
            if self.node_walls[dwn1,nw]==-1:
                self.node_walls[dwn1,nw] = wallSeg
                self.node_walls[dwn1,nw+1] = -1
                break
        for nw in range(self.maxNbs):
            if self.node_walls[dwn2,nw]==-1:
                self.node_walls[dwn2,nw] = wallSeg
                self.node_walls[dwn2,nw+1] = -1
                break
        daughterwalls0.append(wallSeg)
        daughterwalls1.append(wallSeg)

        # create new division wall
        #divWall = self.nextWallIdx
        #self.nextWallIdx += 1
        #self.wall_nodes[divWall] = [divWallNodes[0], divWallNodes[1]]
        #self.wall_cells[divWall] = [cidx, newCell]
        #self.wall_stiff[divWall] = 1.0

        # add div wall to its nodes
        #for n in divWallNodes:
        #    for nw in range(self.maxNbs):
        #        if self.node_walls[n,nw]==-1:
        #            self.node_walls[n,nw] = divWall
        #            self.node_walls[n,nw+1] = -1
        #            break

        # sort walls by which side of division axis they lie
        #daughterwalls0=[divWall]
        #daughterwalls1=[divWall]
        for w in self.cell_walls[cidx,:]:
            if w==-1:
                break
            # This is annoying - pyopencl float4 type does not 
            # initialise numpy.array (like a tuple e.g.), so have 
            # to put it into a list then array
            np0=self.nodes[self.wall_nodes[w,0]] 
            p0 = numpy.array([np0[i] for i in range(2)])
            np1=self.nodes[self.wall_nodes[w,1]] 
            p1 = numpy.array([np1[i] for i in range(2)])
            p = 0.5*(p0+p1) # wall centre
            #print "wall centre " + str(w) + " = "
            #print p
            side=self.whichSideWall(wallSeg, p)
            #print "side = " + str(side)
            if side<0.0:
                # old cell
                daughterwalls0.append(w)
            #elif side<0.0:
            else:
                # new cell
                daughterwalls1.append(w)
                if self.wall_cells[w,0]==cidx:
                    self.wall_cells[w,0] = newCell
                else:
                    self.wall_cells[w,1] = newCell

        # Set walls of two daughter cells
        self.cell_walls[cidx,0:len(daughterwalls0)] = daughterwalls0
        self.cell_walls[cidx,len(daughterwalls0)] = -1 #end of list
        self.cell_walls[newCell,0:len(daughterwalls1)] = daughterwalls1
        self.cell_walls[newCell,len(daughterwalls1)] = -1 #end of list

        # Set turgor of new cell
        self.cell_turgor[newCell] = self.P

    def synchDevData(self):
        # Copy to device (only need to do this max once per step)
        numCells = self.nextCellIdx
        numNodes = self.nextNodeIdx
        numWalls = self.nextWallIdx
        cl.enqueue_copy(self.queue, self.nodes_dev.data, self.nodes[0:numNodes])
        cl.enqueue_copy(self.queue, self.disp_dev.data, self.disp[0:numNodes])
        cl.enqueue_copy(self.queue, self.wall_nodes_dev.data, self.wall_nodes[0:numWalls,:])
        cl.enqueue_copy(self.queue, self.wall_cells_dev.data, self.wall_cells[0:numWalls,:])
        cl.enqueue_copy(self.queue, self.wall_stiff_dev.data, self.wall_stiff[0:numWalls])
        cl.enqueue_copy(self.queue, self.node_walls_dev.data, self.node_walls[0:numNodes,:])
        cl.enqueue_copy(self.queue, self.cell_walls_dev.data, self.cell_walls[0:numCells,:])
        cl.enqueue_copy(self.queue, self.cell_turgor_dev.data, self.cell_turgor[0:numCells])
        self.queue.finish()

    def energy(self):
        perim = cl_array.max(self.cell_area_dev, self.queue).get()
        A = self.cell_area_dev.get()
        A0 = self.cell_area0_dev.get()
        Adiff = A-A0
        Asq = numpy.dot(Adiff,Adiff)
        energy = self.Karea*Asq + self.Kperim*perim
        return energy

    def gradientDescent(self, step):
        # Update data on cell shapes etc.
        self.updateGeometry()
        #print "Energy = %f"%(self.energy())

        # Compute dE/dx
        self.FEMKernels.energyDeriv(self.queue, (self.nextNodeIdx,), None,
                            self.nodes_dev.data,
                            self.node_walls_dev.data,
                            self.wall_nodes_dev.data,
                            self.wall_cells_dev.data,
                            self.cell_area_dev.data,
                            self.cell_area0_dev.data,
                            self.cell_perim_dev.data,
                            self.energy_deriv_dev.data,
                            numpy.float32(self.Karea), numpy.float32(self.Kperim),
                            self.maxNbs).wait()
        # Add step * -dE/dx to node pos
        self.vsubkx(self.nodes_dev, numpy.float32(step), self.nodes_dev, self.energy_deriv_dev)

    def Ax(self, disp, result, dt):
        numNodes = (self.nextNodeIdx,)
        # force=Ax, x=disp
        self.FEMKernels.Ax(self.queue, numNodes, None, \
                self.nodes_dev.data, disp.data,\
                self.node_walls_dev.data, self.wall_nodes_dev.data, \
                self.wall_stiff_dev.data, self.wall_length_dev.data, result.data, \
                self.S, self.I, numpy.float32(dt), \
                self.maxNbs, self.talpha).wait()
	#su = cl_array.sum(disp).get()
	#print "su = "
	#print su
	#self.vaddfloat4(result, su, result)
	
    def CGSSolve(self, dt):
        # Solve Ku=f (Ax=b)
        # u = displacements, in disp_dev device array

        # There must be a way to do this using built in pyopencl - what is it?!
        self.vclear(self.rhs_dev)

        numNodes = (self.nextNodeIdx,)

        # put force in rhs (b)
        self.FEMKernels.rhs(self.queue, numNodes, None, \
                self.nodes_dev.data, self.node_walls_dev.data, \
                self.wall_cells_dev.data, self.wall_nodes_dev.data, \
                self.cell_turgor_dev.data, self.rhs_dev.data, 
		self.disp_dev.data, self.talpha, self.maxNbs).wait()
	# Reset disp
        self.vclear(self.disp_dev)


        # force=Ax, x=disp
        #self.FEMKernels.Ax(self.queue, numNodes, None, \
        #        self.nodes_dev.data, self.disp_dev.data,\
        #        self.node_walls_dev.data, self.wall_nodes_dev.data, \
        #        self.wall_stiff_dev.data, self.wall_length_dev.data, self.force_dev.data, \
        #        self.S, self.I, numpy.float32(dt), \
        #        self.maxNbs, self.talpha).wait()
        self.Ax(self.disp_dev, self.force_dev, dt)
	# res = b-Ax
        self.vsub(self.res_dev, self.rhs_dev, self.force_dev)
        #print "force (Ax) = "
                #print self.force_dev.get()[0:self.nextNodeIdx]
        #print "rhs (b) = "
        #print self.rhs_dev.get()
        #print "res_dev = "
        ##print self.res_dev.get()

        # p = res
        cl.enqueue_copy(self.queue, self.p_dev.data, self.res_dev.data)
        # rsold = l2norm(res)
        rsold = self.vdot(self.res_dev,self.res_dev).get()
        #print "Initial l2norm(res) = " + str(rsold)

        # iterate
        # max iters = matrix dimension = 2 (dofs) * num nodes
        nnodes = self.nextNodeIdx
        for iter in range(nnodes*2):
            # Ap
            #self.FEMKernels.Ax(self.queue, numNodes, None, \
            #    self.nodes_dev.data, self.p_dev.data,\
            #    self.node_walls_dev.data, self.wall_nodes_dev.data, \
            #    self.wall_stiff_dev.data, self.wall_length_dev.data, self.Ap_dev.data, \
            #    self.S, self.I, numpy.float32(dt), \
            #    self.maxNbs, self.talpha).wait()
	    self.Ax(self.p_dev, self.Ap_dev, dt)
            #print "self.Ap_dev, iter " + str(iter) + " = "
            #print self.Ap_dev.get()

            # p^TAp
            pAp = self.vdot(self.p_dev, self.Ap_dev).get()
            #print "pAp, iter " + str(iter) + " = " + str(pAp)

            # alpha = rsold/p^TAp
            alpha = numpy.float32(rsold/pAp)
            #print "alpha, iter " + str(iter) + " = " + str(alpha)

            # x = x + alpha*p, x=self.disp
            self.vaddkx(self.disp_dev, alpha, self.disp_dev, self.p_dev)

            # res = res - alpha*Ap
            self.vsubkx(self.res_dev, alpha, self.res_dev, self.Ap_dev)
            #print "self.res_dev, iter " + str(iter) + " = "
            #print self.res_dev.get()
            #print "self.p_dev, iter " + str(iter) + " = "
            #print self.p_dev.get()

            # rsnew = l2norm(res)
            rsnew = self.vdot(self.res_dev, self.res_dev).get()

            # Test for convergence
            if math.sqrt(rsnew/nnodes)<self.cgsTol:
                break

            # p = res + rsnew/rsold *p
            self.vaddkx(self.p_dev, numpy.float32(rsnew/rsold), self.res_dev, self.p_dev)
            rsold = rsnew
            print "Residual at iteration " + str(iter) + " = " + str(rsnew)

        #print self.disp_dev.get()
        """self.FEMKernels.Ax(self.queue, numNodes, None, \
                        self.nodes_dev.data, self.disp_dev.data,\
                self.node_walls_dev.data, self.wall_nodes_dev.data, \
                self.wall_stiff_dev.data, self.force_dev.data, \
                self.maxNbs, self.talpha).wait()
                """
        #print "final self.force = "
        #print self.force_dev.get()
        #print "rhs = "
        #print self.rhs_dev.get()
                #print numNodes
        #self.vsub(self.res_dev, self.rhs_dev, self.force_dev)
        #print "Final residual calc = "
        #print self.res_dev.get()

        print "CGS residual = %g, rms res = %g, num nodes = %d, ncells = %d, iterations = %d"%(rsnew,math.sqrt(rsnew/nnodes),nnodes,self.nextCellIdx,iter)

    def updateGeometry(self):
        # Compute wall lengths
        self.FEMKernels.computeWallLength(self.queue, (self.nextWallIdx,), None, \
                self.wall_nodes_dev.data, self.nodes_dev.data,  
                self.wall_length_dev.data).wait()

        #l = self.wall_length_dev.get()
        #print "len min = %f"%numpy.min(l[0:model.nextWallIdx])
        #print "len max = %f"%numpy.max(l[0:model.nextWallIdx])
        #print l[0:self.nextWallIdx]
        self.FEMKernels.computeCellGeometry(self.queue, (self.nextCellIdx,), None,\
                        self.cell_walls_dev.data,\
                        self.wall_nodes_dev.data, \
                        self.wall_cells_dev.data, \
                        self.nodes_dev.data, \
                        self.cell_area_dev.data, \
                        self.cell_centre_dev.data, \
                        self.cell_perim_dev.data, \
                        numpy.int32(self.maxCellWalls)).wait()
        #print "Cell areas = "
        self.cell_area = self.cell_area_dev.get()
        #print cellArea
        #print "Cell centres = "
        self.cell_centre = self.cell_centre_dev.get()
        #print cellCentre

    def step(self, dt):
        # Update data on cell shapes etc.
        self.updateGeometry()

        # Divide any cells with area>200
        # ---
        maxArea = cl_array.max(self.cell_area_dev, self.queue).get()
        print "max area = " + str(maxArea)
        if maxArea>400.0:
            self.FEMKernels.computeCellAxis(self.queue, (self.nextCellIdx,), None,\
                    self.cell_walls_dev.data, \
                    self.wall_nodes_dev.data, \
                    self.nodes_dev.data, \
                    self.cell_centre_dev.data, \
                    self.cell_axis_dev.data, \
                    numpy.int32(self.maxCellWalls)).wait()
            self.cell_axis = self.cell_axis_dev.get()
            self.nodes = self.nodes_dev.get()
            self.disp = self.disp_dev.get()
            #print "Cell axes= "
            #print cellAxis 
            nCells = self.nextCellIdx
            for c in range(nCells):
                #print "Cell %d area = %f"%(c,self.cell_area[c])
                if self.cell_area[c]>400.0:
                    ctr= self.cell_centre[c] 
                    ctra = numpy.array([ctr[i] for i in range(2)])
                    axis= self.cell_axis[c] 
                    axisa = numpy.array([axis[i] for i in range(2)])
                    self.divide_cell(c, 0, 0, ctra, axisa) 
            self.synchDevData()
            self.updateGeometry()
            #printAll()
        #---


        # Solve for displacements
        self.CGSSolve(dt)
        # Add self.displacement solution *dt to node positions
        self.vaddkx(self.nodes_dev, numpy.float32(1.0), self.nodes_dev, self.disp_dev)
        
        # force=Ax, x=disp
        #model.FEMKernels.Ax(model.queue, (model.nextNodeIdx,), None, \
        #        model.nodes_dev.data, model.disp_dev.data,\
        #        model.node_walls_dev.data, model.wall_nodes_dev.data, \
        #        model.wall_stiff_dev.data, model.wall_length_dev.data, model.rhs_dev.data, \
        #        model.S, model.I, numpy.float32(dt), \
        #        model.maxNbs, model.talpha).wait()

        # check result
        #print "Node pos "
                #self.nodes = self.nodes_dev.get()
        #print self.nodes

        if self.useVBO:
            cl.enqueue_acquire_gl_objects(self.queue, [self.coords_dev, self.idx_dev])
            cl.enqueue_copy(self.queue, self.coords_dev, self.nodes_dev.data, byte_count=4*4*self.nextNodeIdx)
            cl.enqueue_copy(self.queue, self.idx_dev, self.wall_nodes_dev.data, byte_count=2*4*self.nextWallIdx)
            cl.enqueue_release_gl_objects(self.queue, [self.coords_dev, self.idx_dev])
            self.queue.finish()

    def matrixTest(self):
        x_dev = cl_array.zeros(self.queue, (self.nextNodeIdx,), cl_array.vec.float4)
        Ax_dev = cl_array.zeros(self.queue, (self.nextNodeIdx,), cl_array.vec.float4)
        opstring = ''
        for i in range(self.nextNodeIdx):
            x = numpy.zeros((self.nextNodeIdx,), cl_array.vec.float4)
            for j in range(3):
                if j>0:
                    x[i][j-1]=0.0
                x[i][j]=1.0
                x_dev.set(x)
                self.FEMKernels.Ax(self.queue, (model.nextNodeIdx,), None, \
                    self.nodes_dev.data, x_dev.data,\
                    self.node_walls_dev.data, self.wall_nodes_dev.data, \
                    self.wall_stiff_dev.data, self.wall_length_dev.data, Ax_dev.data, \
                    self.S, self.I, numpy.float32(1), \
                    self.maxNbs, self.talpha).wait()
                #self.calculate_Ax(Ax_dev, x_dev)
                Ax = Ax_dev.get()
                for ii in range(self.nextNodeIdx):
                    for jj in range(3):
                        opstring += '%e'%(Ax[ii][jj])
                        if ii!=self.nextNodeIdx-1 or jj!=2:
                            opstring = opstring + '\t'
                opstring = opstring + '\n'
        print "Ax"
        print opstring
        open('matrix.mat', 'w').write(opstring)


def printAll():
    print "nodes = "
    print model.nodes[0:model.nextNodeIdx]
    print "wall_nodes = "
    print model.wall_nodes[0:model.nextWallIdx]
    print "wall_cells = "
    print model.wall_cells[0:model.nextWallIdx]
    print "node_walls = "
    print model.node_walls[0:model.nextNodeIdx]
    print "cell_walls = "
    print model.cell_walls[0:model.nextCellIdx]

def display():
    glEnable(GL_POINT_SMOOTH)
    glEnable(GL_LINE_SMOOTH)
    glEnable(GL_POLYGON_SMOOTH)
    glEnable(GL_BLEND)
    glEnable(GL_CULL_FACE)
    glCullFace(GL_BACK)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glClearColor(1, 1, 1, 0.7)
    glClear(GL_COLOR_BUFFER_BIT)
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glScalef(0.001,0.001,0.001)
    #glTranslatef(-50,-50,0)
    glPointSize(1.0)
    glLineWidth(2.0)
    glColor3f(0.0,0.0,0.0)

    if model.useVBO:
        glColor3f(1,1,1)
        glBegin(GL_POINTS)
        glVertex3fv([0,0,0])
        glEnd()
        glColor3f(0.0,0.0,0.0)
        glEnableClientState(GL_VERTEX_ARRAY)
        glVertexPointer(3, GL_FLOAT, 4*4, None)
        glBindBuffer(GL_ARRAY_BUFFER, model.vbo)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model.ibo)
        #glDrawArrays(GL_POINTS, 0, model.nextNodeIdx)
        glDrawElements(GL_LINES, model.nextWallIdx*2, GL_UNSIGNED_INT, None)
        glFlush()
        glDisableClientState(GL_VERTEX_ARRAY)
    else:
        glColor3f(0.5,0.5,0.5)
        glBegin(GL_TRIANGLES)
        for i in range(model.nextWallIdx):
            idx = model.wall_nodes[i,:]
            c1 = model.wall_cells[i,0]
            c2 = model.wall_cells[i,1]
            cp1 = list(model.cell_centre[c1])
            cp2 = list(model.cell_centre[c2])
            p1 = list(model.nodes[idx[0]])
            p2 = list(model.nodes[idx[1]])
            glVertex2fv(cp1[0:2])
            glVertex2fv(p2[0:2])
            glVertex2fv(p1[0:2])
            #glVertex2fv(cp2[0:2])
            #glVertex2fv(p1[0:2])
            #glVertex2fv(p2[0:2])
        glEnd()
        glColor3f(0.0,0.0,0.0)
        glLineWidth(2.0)
        glBegin(GL_LINES)
        for i in range(model.nextWallIdx):
            idx = model.wall_nodes[i,:]
            p1 = list(model.nodes[idx[0]])
            p2 = list(model.nodes[idx[1]])
            glVertex2fv(p1[0:2])
            glVertex2fv(p2[0:2])
        #p1 = [0,0]
        #p2 = [1,0]
        glEnd()
        glBegin(GL_POINTS)
        for i in range(model.nextNodeIdx):
            p = list(model.nodes[i])
            glVertex2fv(p[0:2])
        glColor3f(1,0,0)
        for i in range(model.nextCellIdx):
            p = list(model.cell_centre[i])
            glVertex2fv(p[0:2])
        #p1 = [0,0]
        #p2 = [1,0]
        glEnd()
        
        #rhs = model.rhs_dev.get()
        #glColor3f(1.0,0.0,0.0)
        #glLineWidth(1.0)
        #glBegin(GL_LINES)
        #for i in range(model.nextNodeIdx):
        #    p1 = list(model.nodes[i])
        #    f = list(rhs[i])
        #    print "f = " + str(f)
        #    glVertex2fv(p1[0:2])
        #    glVertex2f(p1[0]+2e4*f[0], p1[1]+2e4*f[1])
        #glEnd()


    glFlush()
    glPopMatrix()
    #glFinish()

def reshape(w, h):
    glViewport(0, 0, w, h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glMatrixMode(GL_MODELVIEW)

def animate():
    model.step(1800)
    #if model.nextCellIdx>1:
    #    model.matrixTest()
    '''print "Step"
    print "----"
    model.updateGeometry()
    # Set target area to 1.2*area
    model.cell_area0_dev.set( model.cell_area_dev.get() + 1e3 )
    for i in range(10):
        model.gradientDescent(1e-5)'''
    model.nodes = model.nodes_dev.get()
    model.cell_centre = model.cell_centre_dev.get()
    glutPostRedisplay()

# The model object
#model = None

if __name__ == '__main__':
    import sys
    glutInit(sys.argv)
    glutInitWindowSize(800, 800)
    glutInitWindowPosition(0, 0)
    glutCreateWindow('fem test')
    glutDisplayFunc(display)
    glutReshapeFunc(reshape)
    glutIdleFunc(animate)
    # S=1.0, I=0.0, mu=100.0, P=0.00001, 
    model = CLPlant(1e-3, 1e-4, useVBO=True, I=0)
    #model.divide_cell(0, 2,3, numpy.array([0.5,0.5]), numpy.array([1.0,0.0]))
    print "nodes = "
    print model.nodes
    print "wall_nodes = "
    print model.wall_nodes
    print "wall_cells = "
    print model.wall_cells
    print "node_walls = "
    print model.node_walls
    print "cell_walls = "
    print model.cell_walls

    import timeit
    s = '''
from __main__ import CLPlant
model = CLPlant(0.5, 1e-6)
        '''
#        t = timeit.Timer(stmt='model.step(0.2)', setup=s)
#        try:
#                print "Time = "
#                print t.timeit(5000)/5000
#        except:
#                print t.print_exc()

    axis = numpy.array([1.0,0.0])
    for i in range(0):
        numCells = model.nextCellIdx
        for c in range(numCells):
            Ltot=0.0
            avg = numpy.array([0,0,0])
            wi=0
            while model.cell_walls[c,wi]!=-1:
                wall = model.cell_walls[c,wi]
                p0 = numpy.array([model.nodes[model.wall_nodes[wall,0]][i] for i in range(2)])
                p1 = numpy.array([model.nodes[model.wall_nodes[wall,1]][i] for i in range(2)])
                L = numpy.linalg.norm(p0-p1)
                Ltot += L
                avg = avg + 0.5*(p0+p1)*L 
                wi += 1
            avg = avg / Ltot #numpy.float32(wi)

            model.divide_cell(c, 0,0, avg[0:2]+ numpy.random.uniform(-0.1,0.1,2) , axis)
            a0 = axis[0]
            axis[0] = axis[1]
            axis[1] = a0
            print "numCells = " + str(model.nextCellIdx)

    model.synchDevData()

    print "numNodes = " + str(model.nextNodeIdx)

    #for i in range(1):
    #        model.step(1800)
    #model.nodes = model.nodes_dev.get()
    #print "nodes = "
    #print model.nodes[0:model.nextNodeIdx]
    #print "wall_nodes = "
    #print model.wall_nodes[0:model.nextWallIdx]
    #print "wall_cells = "
    #print model.wall_cells[0:model.nextWallIdx]
    #print "node_walls = "
    #print model.node_walls[0:model.nextNodeIdx]
    #print "cell_walls = "
    #print model.cell_walls[0:model.nextWallIdx]
    #print "cell_nodes = "
    #print model.cell_nodes[0:model.nextCellIdx]

    glutMainLoop()

    #model.updateGeometry()
    #model.CGSSolve(1.0)
    #model.nodes = model.nodes_dev.get()
    #print "nodes = "
    #print model.nodes
    #model.disp = model.disp_dev.get()
    #print "disp = "
    #print model.disp
    #model.rhs = model.rhs_dev.get()
    #print "rhs = "
    #print model.rhs

