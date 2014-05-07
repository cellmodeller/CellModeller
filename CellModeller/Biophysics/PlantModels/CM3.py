#
# CellModeller3
#
# core module
#
# Jonathan Mackenzie
# Lionel Dupuy  
#
# Nov 2007
#

# CellModeller4 interface
# Tim Rudge
# Jan 2012

# Python 
import random 
import copy
import math
# site_packages
import numpy
from scipy import weave, sparse, linalg
import scipy.sparse.linalg.dsolve as linsolve
import xml.dom.minidom as minidom
import sys

class CM3:
    """
    Wraps all that is needed for a simulation,
    and defines interface for interacting with it.
    """
    def __init__(self, sim):
        # Simulator object
        self.sim = sim

        # object counters
        self.numnodes=0 
        self.numwalls=0
        self.numcells=0
        self.celltotal=0
        self.numcelltypes=0
        
        self.NUMDOF=3 # number of degrees of freedom of each node
        # Finite Element Model constants (from elastic.cpp)
        self.FEM_S=1.0e-11 # wall cross section, m^2
        self.FEM_E=6.0e12 # viscosity coefficient, GPa.s (Giga Pascal secs - 1 Pa-1N/m^2)
        self.FEM_nu=0.3 # Poisson modulus, no dimensions (how a solid deforms in the y/z direction if you stretch it in the x direction e.g. an elastic)
        #self.FEM_G=FEM_E/(2.0*(1.0+FEM_nu)) # shear modulus - not used in 2D
        self.FEM_I=0 #5.0e-21 # moment of inertia or quadratic moment, m^4
        #self.s_scale = 1e-6 # spatial scale, micrometer
        #self.t_scale = 3600.0 # time scale = 1 hour
    
        # numpy array-based data structures
        self.nodedisparray=numpy.empty((0,1)) # node displacement column vector - NUMDOF elements per node stacked vertically
        self.nodeforcearray=numpy.empty((0,1)) # node forces column vector - NUMDOF elements per node stacked vertically
        
        # linked-list data structures which map onto numpy data
        self.nodes=[] # list of Node objects
        self.walls=[] # list of Wall objects
        #self.cells.values()=[] # list of Cell objects
        self.cells = {} # map from cell_id to cell
        
        self.celltypes=[] # list of CellType objects
        ct = CellType(self)
        self.nodestofix=[]
        
        #self.create_one_cell()
        #self.create_cell_grid()
        
        self.stiffness = numpy.zeros((self.numnodes*self.NUMDOF,self.numnodes*self.NUMDOF))
    
        def reset(self):
        self.__init__(self.sim)

    def setRegulator(self, reg):
        self.reg = reg

    def hasNeighbours(self):
        return True

    def divide(self, parentState, d1State, d2State, *arg, **kw):
        cell = self.cells[parentState.id]
        
        #cell.calc_principle_axes()
        axis=cell.pa2
        d1_id = d1State.id
        d2_id = d2State.id
        self.divide_cell(cell, d1_id, d2_id, cell.get_centre(), axis)
        for cell in self.cells.values():
            cell.order_nodes()

        # Initialise new cell state
        self.initCellState(d1State)
        # Update neighbourhoods
        for n in d1State.nbs:
            self.updateNeighbourhood(self.sim.cellStates[n])
        # Initialise new cell state
        self.initCellState(d2State)
        # Update neighbourhoods
        for n in d2State.nbs:
            self.updateNeighbourhood(self.sim.cellStates[n])

    
    def divide_cell(self, cell, d1id, d2id, centre, axis):
        """
        Divide cell along axis passing through centre
        """
        # find wall intersections
        intersections=[]
        for w in cell.walls:
            p=w.intersect(centre,axis)
            if p!=None:
                intersections.append((w,p))
                
        if len(intersections)!=2:
            raise Exception,"Cell division axis does not divide 2 walls"
        
        # remove intersected walls, add new split walls, and refs 
        newnodes=[]
        for w,p in intersections:
            neighbourcell=w.get_neighbour(cell)
            self.remove_wall(w)
            
            # add new node and refs
            newnode=self.create_node()
            newnode.set_pos(p)
            newnodes.append(newnode)
            if neighbourcell:
                neighbourcell.nodes.append(newnode)
                newnode.cells.append(neighbourcell)
                
            # add new walls
            for n in (w.node0,w.node1):
                newwall=self.create_wall(n,newnode) 
                cell.walls.append(newwall)
                # add refs to neighbour, if it exists
                if neighbourcell:
                    neighbourcell.walls.append(newwall)
                    newwall.append_cell(neighbourcell)
                    #neighbourcell.order_nodes()

                        
        # create new division wall
        newwall=self.create_wall(newnodes[0],newnodes[1])
        
        # sort walls by which side of division axis they lie
        daughterwalls0=[newwall]
        daughterwalls1=[newwall]
        for w in cell.walls:
            #print( w,newwall.which_side(w.node0.get_pos()),newwall.which_side(w.node1.get_pos()) )
            side=newwall.which_side(w.get_centre())
            if side>0.0:
                daughterwalls0.append(w)
            #elif side<0.0:
            else:
                daughterwalls1.append(w)
                                    
        # remove mother cell, and add two new daughter cells
        mothertype=cell.type
        self.remove_cell(cell)
        daughter=self.create_cell(daughterwalls0, d1id)
        daughter.type=mothertype
        daughter=self.create_cell(daughterwalls1, d2id)
        daughter.type=mothertype

    def step(self, dt):
        self.tick_fem(dt*25.0)
        self.updateCellStates()
        
    def dp(self):
        # debug print
        for w in self.walls:
            print w
        for n in self.nodes:
            print n

        def updateCellStates(self):
        nbmap = self.getNeighbours()
        for (cid,cellState) in self.sim.cellStates.items():
            cell = self.cells[cid]
            cellState.pos = cell.get_centre()
            cellState.volume = cell.get_area()
            cellState.nbs = nbmap[cid]
            cellState.wallp = []
            for w in cell.walls:
                cellState.wallp.append( (w.node0.get_pos(), w.node1.get_pos()) )
            cellState.nodep = []
            for n in cell.nodes:
                cellState.nodep.append( n.get_pos() )
            # Growth rate?
            cell.growthrate = cellState.growthRate

    def updateNeighbourhood(self, cellState):
        cell = self.cells[cellState.id]
        cellState.nbs = self.getCellNeighbours(cell)

    def initCellState(self, cellState):
        cell = self.cells[cellState.id]
        cellState.pos = cell.get_centre()
        cellState.volume = cell.get_area()
        cellState.nbs = self.getCellNeighbours(cell)
        cellState.startVol = cellState.volume
        cellState.wallp = []
        for w in cell.walls:
          cellState.wallp.append( (w.node0.get_pos(), w.node1.get_pos()) )
        cellState.nodep = []
        for n in cell.nodes:
            cellState.nodep.append( n.get_pos() )
        # Growth rate?
            
    def addCell(self, cellState):
        # add an N sided cell
        N=5
        #o=math.pi/4.0
        o=0.0
        #s=100.0/0.7071
        s=10.0
        for i in range(N):
            node=self.create_node()
            node.set_pos((s*math.cos(i*2*math.pi/N+o),s*math.sin(i*2*math.pi/N+o)))
            
        for i in range(N):
            wall=self.create_wall(self.nodes[i],self.nodes[(i+1)%N])
        
        cell=self.create_cell(self.walls, cellState.id)
        self.initCellState(cellState)
                    
    def fix(self, tofix=True):
        # fix or unfix nodes of picked cells
        self.nodestofix=[]
        for cell in self.pickedcells:
            for node in self.cells.values()[cell].nodes:
                if node not in self.nodestofix:
                    self.nodestofix.append(node)
                                    
        
    def set_cell_types(self, celltypeid):
        # set all selected cells params to that of selected type
        celltype=self.celltypes[celltypeid]
        for cell in self.pickedcells:
            self.cells.values()[cell].type=celltype
            #self.cells.values()[cell].colnorm[0]=celltype.colour.Red()/256.0
            #self.cells.values()[cell].colnorm[1]=celltype.colour.Green()/256.0
            #self.cells.values()[cell].colnorm[2]=celltype.colour.Blue()/256.0
            #self.cells.values()[cell].col=self.cells.values()[cell].colnorm
            #self.cells.values()[cell].turgor=celltype.turgor
            
        
        '''
        ####################################################
        Add and Remove functionality of the model interface.
        In the functions below:
        add and remove functions take care of model id counts
        node, wall and cell constructors take care of adding references
        remove functions take care of removing references
        ####################################################
        '''
                                
    def create_node(self):
        """
        Adds a node to both data structures, and returns a Node instance
        """
        self.nodedisparray=numpy.vstack([self.nodedisparray,numpy.zeros((self.NUMDOF,1))]) # add space to nodedisparray
        self.nodeforcearray=numpy.vstack([self.nodeforcearray,numpy.zeros((self.NUMDOF,1))]) # and to nodeforcearray
        node=Node(self)
        self.nodes.append(node)
        node.id=self.numnodes
        self.numnodes+=1
        return node
    
    def create_wall(self,node0,node1):
        """
        Adds a wall to connect two nodes, and returns a Wall instance
        """
        if node0==node1:
            raise Exception,"Nodes cannot be the same"
        
        wall=Wall(self,node0,node1)
        self.walls.append(wall)
        wall.id=self.numwalls
        self.numwalls+=1

        return wall
    
    def create_cell(self, walls, cid):
        """
        Adds a cell given a list of more than 2 walls, and returns a Cell instance
        """
        if len(walls)<3:
            raise Exception,"Need more than 2 walls to make a cell"
        cell=Cell(self,walls)
        self.cells[cid] = cell
        cell.id= cid #self.celltotal+1
        self.numcells+=1 # this is decremented when cell removed
        self.celltotal+=1 # this never gets decremented -> unique id
                
        return cell
    
    def remove_wall(self,wall):
        """
        Removes a wall and references to it, without removing nodes
        """
        # remove refs to this wall in cells
        if wall.cell0:
            wall.cell0.walls.remove(wall)
        if wall.cell1:
            wall.cell1.walls.remove(wall)
            
        # remove refs to this wall in nodes
        wall.node0.walls.remove(wall)
        wall.node1.walls.remove(wall)
    
        # remove wall from model
        self.walls.remove(wall)
        self.numwalls-=1
        # renumber wall ids
        for w in self.walls:
            if w.id > wall.id:
                w.id-=1
        
    def remove_cell(self,cell):
        """
        Removes a cell and references to it, without removing walls or nodes
        """
        for w in cell.walls:
            if w.cell0==cell:
                w.cell0=None
            if w.cell1==cell:
                w.cell1=None
                
        for n in cell.nodes:
            n.cells.remove(cell)
            
        del self.cells[cell.id]
        self.numcells-=1
        #for c in self.cells.values():
        #    if c.id > cell.id:
        #        c.id-=1
                
        '''
        ###########################################
        ###########################################
        '''
        
    def getNeighbours(self):
        ndict = {}
        for c in self.cells.values():
            nbs = []
            for w in c.walls:
                n = w.get_neighbour(c)
                if n:
                    nbs.append(n.id) 
            ndict[c.id] = nbs
        return ndict        
            
    def getCellNeighbours(self, c):
        nbs = []
        for w in c.walls:
            n = w.get_neighbour(c)
            if n:
                nbs.append(n.id) 
        return nbs

    def tick_fem(self,dt):
        """
        tick Finite Element Model simulation by dt hours
        """
        for cell in self.cells.values():
            cell.calc_principle_axes()        
            
        # initialise stiffness matrix
        self.stiffness=numpy.zeros((self.numnodes*self.NUMDOF,self.numnodes*self.NUMDOF)) 
        #self.stiffness.resize((self.numnodes*self.NUMDOF,self.numnodes*self.NUMDOF)) 
        #self.stiffness[:,:] = 0
        
        # for each wall, add a component to the stiffness matrix
        for wall in self.walls:
            L=wall.get_length()
            
            # Average growth rate of adjacent cells
            g = 0.0
            ncells = 0.0
            if wall.cell0:
                ncells += 1.0
                g += wall.cell0.growthrate
            if wall.cell1:
                ncells += 1.0
                g += wall.cell1.growthrate
            g /= ncells

            Ra = self.FEM_E*self.FEM_S/L/dt
            Rf = self.FEM_E*self.FEM_I/(L**3)/dt
            #Ra = self.FEM_E/g * self.FEM_S/dt #/L/dt
            #Rf = 0.0 #self.FEM_E/g * self.FEM_I/(L**3)/dt
            M =numpy.array(((Ra  , 0   , 0      , -Ra , 0    , 0      ),
                          (0   , Rf*12  , Rf*6*L    , 0   , -Rf*12  , Rf*6*L    ),
                          (0   , Rf*6*L , Rf*4*L**2 , 0   , -Rf*6*L , Rf*2*L**2 ),
                          (-Ra , 0   , 0      , Ra  , 0    , 0      ),
                          (0   , -Rf*12 , -Rf*6*L   , 0   , Rf*12   , -Rf*6*L   ),
                          (0   , Rf*6*L , Rf*2*L**2 , 0   , -Rf*6*L , Rf*4*L**2 )))
#            M =numpy.array(((Ra  , 0   , 0      , -Ra , 0    , 0      ),
#                          (0   , 0  , 0    , 0   , 0  , 0    ),
#                          (0   , 0 , 0 , 0   , 0 , 0 ),
#                          (-Ra , 0   , 0      , Ra  , 0    , 0      ),
#                          (0   , 0 , 0   , 0   , 0   , 0   ),
#                          (0   , 0 , 0 , 0   , 0 , 0 )
#                        ))


            # rotate from local coord frame to global
            #theta=wall.get_angle()

            c=wall.get_cosine() 
            #c=math.cos(theta)
            s=wall.get_sine() 
            #s=math.sin(theta)
            R=numpy.matrix(( (c,-s, 0, 0, 0, 0),
                    (s, c, 0, 0, 0, 0),
                    (0, 0, 1, 0, 0, 0),
                    (0, 0, 0, c,-s, 0),
                    (0, 0, 0, s, c, 0),
                    (0, 0, 0, 0, 0, 1) 
                    ))
            #R=numpy.hstack((R,R))
            #R=numpy.vstack((R,R))
            #print "R = "
            #print R
            
            MR=R*M*R.transpose()
            #print "MR = "
            #print MR
            
            # add component to stiffness matrix
            self.stiffness[wall.node0.id*3:wall.node0.id*3+3,wall.node0.id*3:wall.node0.id*3+3]+=MR[0:3,0:3]
            self.stiffness[wall.node1.id*3:wall.node1.id*3+3,wall.node1.id*3:wall.node1.id*3+3]+=MR[3:6,3:6]
            self.stiffness[wall.node0.id*3:wall.node0.id*3+3,wall.node1.id*3:wall.node1.id*3+3]+=MR[0:3,3:6]
            self.stiffness[wall.node1.id*3:wall.node1.id*3+3,wall.node0.id*3:wall.node0.id*3+3]+=MR[3:6,0:3]

        #print_mat(self.stiffness)
        
        # apply turgor force
        self.nodeforcearray=numpy.zeros(self.nodeforcearray.shape)
        for cell in self.cells.values():
            cell.apply_turgor()            
        #self.cells.values()[0].apply_turgor(0.001)
        
        '''
        solve for displacement (find d in S.d=f)
        '''
        # fix nodes
        if len(self.nodestofix)==0:
            # fixing first N nodes
            N=1
            f=numpy.matrix(self.nodeforcearray) #[3*N:]
            s=numpy.matrix(self.stiffness) #[3*N:,3*N:]
        else:            
            # fix nodes listed in self.nodestofix
            nodeids=[]
            for n in self.nodes:            
                if n not in self.nodestofix:
                    nodeids.append(n.id*self.NUMDOF)
                    nodeids.append(n.id*self.NUMDOF+1)
                    nodeids.append(n.id*self.NUMDOF+2)            
            # trim arrays
            f=self.nodeforcearray[nodeids]
            s=self.stiffness[nodeids,:]
            s=s[:,nodeids]
    
        # solve system
        #print "f = "
        #print f
        # +numpy.matrix(numpy.eye(s.shape[0]))
        (d,res,rnk,sv) =linalg.lstsq(s+numpy.matrix(numpy.eye(s.shape[0])),f, cond=1e-2, overwrite_a=True) #solve(s,f)
        #res = linalg.norm(s*d-f)

        #d = linalg.solve(s.transpose()*s, s.transpose()*f)
        #print "d = "
        #print d
        #print "res = %f"%(res)
        #print "s = "
        #print s
        #if numpy.max(d)>1.0:
    #        sys.exit(0)
        
        # update displacements
        if len(self.nodestofix)==0:
            self.nodedisparray += d #[3*N:]
        else:
            self.nodedisparray[nodeids]+=d

    def tick_fem_sparse(self,dt):
        """
        tick Finite Element Model simulation by dt hours using sparse stiffness matrix
        """
        
        # initialise stiffness matrix
        self.stiffness=sparse.lil_matrix((self.numnodes*self.NUMDOF,self.numnodes*self.NUMDOF))
        #self.stiffness=numpy.zeros((self.numnodes*self.NUMDOF,self.numnodes*self.NUMDOF))
        
        # for each wall, add a component to the stiffness matrix
        for wall in self.walls:
            L=wall.get_length()
                        
            Ra = self.FEM_E*self.FEM_S/L/dt
            Rf = self.FEM_E*self.FEM_I/(L**3)/dt
            M =numpy.array(((Ra  , 0   , 0      , -Ra , 0    , 0      ),
                          (0   , Rf*12  , -Rf*6*L    , 0   , -Rf*12  , -Rf*6*L    ),
                          (0   , -Rf*6*L , Rf*4*L**2 , 0   , Rf*6*L , Rf*2*L**2 ),
                          (-Ra , 0   , 0      , Ra  , 0    , 0      ),
                          (0   , -Rf*12 , Rf*6*L   , 0   , Rf*12   , Rf*6*L   ),
                          (0   , -Rf*6*L , Rf*2*L**2 , 0   , Rf*6*L , Rf*4*L**2 )))
            
            # rotate from local coord frame to global
            theta=wall.get_angle()

            c=wall.get_cosine() #math.cos(theta)
            s=wall.get_sine() #math.sin(theta)
            R=numpy.array(((c,-s,0),
                           (s,c,0),
                           (0,0,1)))
            R=numpy.hstack((R,R))
            R=numpy.vstack((R,R))
            
            MR=R*M*R.transpose()
            
            # add component to stiffness matrix
            self.stiffness[wall.node0.id*3:wall.node0.id*3+3,wall.node0.id*3:wall.node0.id*3+3]+=sparse.lil_matrix(MR[0:3,0:3])
            self.stiffness[wall.node1.id*3:wall.node1.id*3+3,wall.node1.id*3:wall.node1.id*3+3]+=sparse.lil_matrix(MR[3:6,3:6])
            self.stiffness[wall.node0.id*3:wall.node0.id*3+3,wall.node1.id*3:wall.node1.id*3+3]+=sparse.lil_matrix(MR[0:3,3:6])
            self.stiffness[wall.node1.id*3:wall.node1.id*3+3,wall.node0.id*3:wall.node0.id*3+3]+=sparse.lil_matrix(MR[3:6,0:3])

        #print_mat(self.stiffness)
        
        # apply turgor force
        self.nodeforcearray=numpy.zeros(self.nodeforcearray.shape)
        for cell in self.cells.values():
            cell.apply_turgor(0.01)            
        #self.cells.values()[0].apply_turgor(0.001)
        
        # solve for displacement (find d in S.d=f)
        # fixing first node
        f=self.nodeforcearray[3:]
        s=self.stiffness[3:,3:]
        #d=linalg.solve(s,f)
        d=linsolve.spsolve(s.tocsr(),f)
        self.nodedisparray[3:]+=d.reshape(d.shape[0],1)            

    def apply_force_to_node(self,node,force):
        self.nodeforcearray[node.id*self.NUMDOF:node.id*self.NUMDOF+2]+=force    
    
    def load(self,filename,version=2):
        """
        Loads a model from XML file
        version='2' loads & converts geometry from old CellModeller 2 files
        version='3' loads complete model from new CellModeller 3 files
        """
        #create dom
        dom=minidom.parse(filename)
        if version==2:
            # load vertices as nodes
            vertelems=dom.getElementsByTagName("vertex")
            #print vert_elems
            vdic={} # maps cm2 ids in file to cm3 wall instances created here
            for ve in vertelems:
                id=long(ve.getAttribute("_id"))
                pos=eval(ve.getAttribute("pos")) #evaluates to list
                node=self.create_node()
                node.set_pos(pos)
                vdic[id]=node
                
            # load walls
            wallelems=dom.getElementsByTagName("wall")
            wdic={} # maps cm2 ids in file to cm3 wall instances created here
            for we in wallelems:
                id=long(we.getAttribute("_id"))
                vid0=long(we.getAttribute("vert0_id"))
                vid1=long(we.getAttribute("vert1_id"))
                node0=vdic[vid0]
                node1=vdic[vid1]
                wall=self.create_wall(node0,node1)
                wdic[id]=wall
                # for detecting double walls
                neighbour=long(we.getAttribute("neighbour"))
                wall.cm2neighbour=neighbour
                
            # load cells
            cellelems=dom.getElementsByTagName("cell")
            cdic={} # maps cm2 ids in file to cm3 wall instances created here
            for ce in cellelems:
                cid=long(ce.getAttribute("_id"))
                wallelems=ce.getElementsByTagName("wall")
                walls=[]
                for we in wallelems:
                    wid=long(we.getAttribute("_id"))
                    wall=wdic[wid]
                    walls.append(wall)
                cell=self.create_cell(walls)
                cdic[cid]=cell    
                
            # look for all walls having neighbours
            # and eliminate double wall
            for wall in self.walls:
                if wall.cm2neighbour>0: #neighbour exists
                    neighbourwall=wdic[wall.cm2neighbour]
                    wall.cell1=neighbourwall.cell0 # modify refs in wall...
                    wall.cell1.walls.append(wall)# ...and in cell
                    self.remove_wall(neighbourwall) # remove duplicate wall

def rotate_mat(mat,theta):
    # rotates stiffness matrix component coord frame by theta
    c=math.cos(theta)
    s=math.sin(theta)
    R=numpy.array(((c,-s,0),
                   (s,c,0),
                   (0,0,1)))
    matrot=R*mat*R.transpose()
    return matrot

def print_mat(mat):
    # formatted printing of matrix
    x,y=mat.shape
    for i in range(x):
        s='['
        for j in range(y):
            s+="%+1.2f "%mat[i,j]
        s+=']'
        print s

class Node:
    """
    Node class for use in linked list side of cellular data structure
    """
    
    def __init__(self,model):
        self.id=0
        self.walls=[] # list of walls this node is part of
        self.cells=[] # list of cells this node is part of
        
        self.model=model
        
    #def set_pos(self,x,y): # from x and y
        #self.model.nodedisparray[self.id*self.model.NUMDOF]=x
        #self.model.nodedisparray[self.id*self.model.NUMDOF+1]=y

    def set_pos(self,p):
        # set node position
        self.model.nodedisparray[self.id*self.model.NUMDOF]=p[0]
        self.model.nodedisparray[self.id*self.model.NUMDOF+1]=p[1]
        
    def get_pos(self):
        # return position as numpy array
        return self.model.nodedisparray[self.id*self.model.NUMDOF:self.id*self.model.NUMDOF+2]
    
    def set_rot(self,r):
        # set node rotation
        self.model.nodedisparray[self.id*self.model.NUMDOF+2]=r

    def get_rot(self):
        # get node rotation
        r=self.model.nodedisparray[self.id*self.model.NUMDOF+2]
        return (r)
        
class Wall:
    """
    Wall class for use in linked list side of cellular data structure
    """

    def __init__(self,model,node0,node1):
        self.id=0
        self.node0=None # nodes at end of this wall 
        self.node1=None     
        self.cell0=None # cells this wall is part of
        self.cell1=None # one of these can be None for outer cell
        self.model=None
            
        self.model=model
        self.node0=node0
        self.node1=node1
        
        # add refs to this wall
        node0.walls.append(self)
        node1.walls.append(self)    
        
    def append_cell(self,cell):
        # add a cell ref to this wall
        if self.cell0==None:
            self.cell0=cell
        elif self.cell1==None:
            self.cell1=cell
        else:
            raise Exception,"Already two cells to this wall"
        
    def get_neighbour(self,cell):
        # returns cell that neighbours one given through this wall
        # or None if none
        if cell==self.cell0:
            return self.cell1
        elif cell==self.cell1:
            return self.cell0
        else:
            return None
        
    def get_centre(self):
        sum=self.node0.get_pos()+self.node1.get_pos()
        sum/=2.0
        return sum
    
    def get_length(self):
        length=numpy.linalg.norm(self.node1.get_pos()-self.node0.get_pos())
        return length
    
    def get_cosine(self):
        n=self.node1.get_pos()-self.node0.get_pos()
        return float(n[0]/self.get_length())

    def get_sine(self):
        n=self.node1.get_pos()-self.node0.get_pos()
        return float(n[1]/self.get_length())
    
    def get_angle(self):
        # get angle of wall from left to right (ie -pi/2 < angle < pi/2 )
        l=self.get_length()
        if l>0:
            n=self.node1.get_pos()-self.node0.get_pos()
            a=math.atan(n[1]/n[0])
            if n[0]<0:
                a+=math.pi
            return a
        else:
            raise Exception,"Wall length is zero"
        
    def intersect(self,c,v):
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
        
        a=self.node0.get_pos()
        b=self.node1.get_pos()
        d=c+v
        
        den=(b[0]-a[0])*(d[1]-c[1])-(b[1]-a[1])*(d[0]-c[0]) #all 2D positions are numpy array shape(2,1), so second [0] needed
        
        if den==0: # wall and line are parallel
            return None
    
        r=((a[1]-c[1])*(d[0]-c[0])-(a[0]-c[0])*(d[1]-c[1]))/den
    
        tol=0.00 # the smallest tolerance causes small wall segment which create numerical singularity in the physics
        if r>tol and r<(1.0-tol):
            p=a+(b-a)*r
            return p
    
        else:
            return None
        
    def which_side(self,p):
        """
        Determines which side of a wall a point lies by returning 
        the cross product of (p-node0) x (node1 - node0)
        ie +ve one side, -ve the other
        """
        n0=self.node0.get_pos()
        n1=self.node1.get_pos()
        
        #x=(p[1]-n0[1])*(n1[0]-n0[0])-(p[0]-n0[0])*(n1[1]-n0[1])
        #return x[0]
        return numpy.cross((p-n0).transpose(),(n1-n0).transpose())[0]

    def get_perp(self,cell):
        # get a unit vector perpendicular to wall, facing outwards from centre of cell
        vp=self.node1.get_pos()-self.node0.get_pos()
        vp/=self.get_length() #normalise
        vpr=numpy.array((-vp[1],vp[0])) #rotate 90 anticlockwise
        p=self.node0.get_pos()+vpr
        # if on same side of wall as centre of cell, flip
        if self.which_side(p) * self.which_side(cell.get_centre()) > 0.0:
            vpr*=-1
        return vpr
    
class Cell:
    """
    Cell class for use in linked list side of cellular data structure
    """

    def __init__(self,model,walls):
        self.id=-1 # unique id and index into model.cells[]
        self.walls=[] # list of walls that makes this cell
        self.nodes=[] # list of nodes that makes this cell
        
        self.model=model
        self.walls=copy.copy(walls)
        self.order_nodes() #creates nodes list from walls
        
        # principle axes
        self.pa1=numpy.zeros((2,1))
        self.pa2=numpy.zeros((2,1))
        
        # add refs to this cell
        # in walls
        for w in self.walls:
            if w.cell0==None:
                w.cell0=self
            elif w.cell1==None:
                w.cell1=self
            else:
                raise Exception,"Wall already belongs to 2 cells"
            
        # in nodes
        for n in self.nodes:
            n.cells.append(self)
            
                
        # cell type
        self.type=self.model.celltypes[0]
        
        # colour
        #self.col=self.type.get_colour_for_gl()

        # growth rate (relative)
        self.growthrate = 1.0
            
    def getId(self):
        return self.id
    def setId(self, id):
        self.id = id
        
    def pos(self):
        return self.get_centre().transpose().tolist()[0]

    def setGrowthRate(self, r):
        self.growthrate = r
        
    def print_to_string(self):
        text='Cell: '+str(self.id)+'\n'
        text+='Walls: '
        for w in self.walls:
            text+=str(w.id)+' '
        text+='\nNodes: '
        for n in self.nodes:
            text+=str(n.id)+' '
            
        return text
        
    def order_nodes(self):
        # walk through walls to create ordered set of nodes
        self.nodes=[]
        currentwall=self.walls[0]
        firstnode=currentwall.node0
        self.nodes.append(firstnode)
        currentnode=currentwall.node1
        while currentnode!=firstnode:
            self.nodes.append(currentnode)
            found=False
            for trialwall in self.walls:
                if trialwall!=currentwall:
                    if trialwall.node0==currentnode:
                        currentnode=trialwall.node1
                        currentwall=trialwall
                        found=True
                        break
                    elif trialwall.node1==currentnode:
                        currentnode=trialwall.node0
                        currentwall=trialwall
                        found=True
                        break
            if not found:
                raise Exception,"Cell nodes not in loop"
            
    def get_centre(self):
        # get average of node positions
        # we need this because nodes are not ordered - a pain
        avg=numpy.zeros((2,1))
        for n in self.nodes:
            avg+=n.get_pos()
        avg/=len(self.nodes)

        # calc CofM of cell as
        # area weighted average of node pos
        pos = numpy.zeros((2,1))
        atot=0
        numnodes=len(self.nodes)
        for i in range(numnodes):
            n0=self.nodes[i]
            n1=self.nodes[(i+1)%numnodes]
            p0=n0.get_pos()-avg
            p1=n1.get_pos()-avg
            a = abs(p0[0]*p1[1]-p1[0]*p0[1]) # abs to make sure of ordering
            pos += a*(p0+p1)
            atot += a
            
        return avg + pos/(atot*6.0)

    #vec3 CCell::calcCentreW()
#{
#// centre according to walls (wall centres, weighted by wall lengths)

    #vec3 p(0,0,0);
    #vec3 q;
    #float totlength=0;
    #std::vector<CWall*>::iterator wit;
    #for(wit=walls.begin(); wit!=walls.end(); wit++)
    #{
        #CWall* w=*wit;
        #q=w->vert0->p + w->vert1->p;
        #float l=w->getLength();
        #totlength+=l;
        #q*=0.5*l;
        #p+=q;
    #}

    #//if(walls.size()>0)
    #//p/=((float)walls.size())*totlength);
    #p/=totlength;
    #centre=p;
    #return p;
#}
    def apply_turgor(self,pressure=None):
        # add forces to nodes perpendicular to wall and proportional to wall length
        # if no pressure given, use cell's type value of turgor
        if not pressure:
            pressure=self.type.turgor
        for wall in self.walls:
            vp=wall.get_perp(self)
            force=vp*pressure*wall.get_length()/2.0
            self.model.apply_force_to_node(wall.node0,force)
            self.model.apply_force_to_node(wall.node1,force)
            
    def calc_principle_axes(self):
        # calculate principle axes
        # see www.cs.princeton.edu/courses/archive/fall03/cs597D/lectures/rigid_registration.pdf
        
        # (1) calc covariance matrix of each wall centre relative to centre of cell, weighted by wall length
        #c=self.calcCentreW()
        cov=numpy.array([[0.0,0.0],[0.0,0.0]])
        c=self.get_centre()
        
        for v in self.nodes:
            q=v.get_pos()-c
            cov[0,0]+=q[0]*q[0]
            cov[0,1]+=q[0]*q[1]
            cov[1,0]+=q[1]*q[0]
            cov[1,1]+=q[1]*q[1]            
        cov/=len(self.nodes)

        #for w in self.walls:
            #q=w.vert0.p+w.vert1.p
            #q*=0.5
            #q-=c
            #cov[0,0]+=q[0]*q[0]*w.getLength()
            #cov[0,1]+=q[0]*q[1]*w.getLength()
            #cov[1,0]+=q[1]*q[0]*w.getLength()
            #cov[1,1]+=q[1]*q[1]*w.getLength()
        #cov/=len(self.walls)
        
        # (2) get eigenvectors and eigenvalues
        (w,vr)=numpy.linalg.eig(cov)
        w/=numpy.sum(w) # normalise
        w=abs(w) 
        
        # (3) calc eigenvectors from centre scaled by eigenvalues and make pa1 the largest 
        if w[0]>w[1]:
            self.pa1=vr[:,0].reshape((2,1))#*w[0]
            self.pa2=vr[:,1].reshape((2,1))#*w[1]
        else:
            self.pa2=vr[:,0].reshape((2,1))#*w[0]
            self.pa1=vr[:,1].reshape((2,1))#*w[1]
        
    def volume(self):
        return self.get_area()
        
    def get_area(self):
        # calc area of cell
        # see http://mathworld.wolfram.com/PolygonArea.html
        
        a=0
        numnodes=len(self.nodes)
        for i in range(numnodes):
            n0=self.nodes[i]
            n1=self.nodes[(i+1)%numnodes]
            p0=n0.get_pos()
            p1=n1.get_pos()
            a+=p0[0]*p1[1]-p1[0]*p0[1]
            
        return abs(a[0]*0.5)


        
class CellType:
    '''
    Cell type data structure 
    '''
    def __init__(self,model):
        # create a cell type, and add it to the model
        self.model=model
        self.id=model.numcelltypes
        model.numcelltypes+=1
        model.celltypes.append(self)
        
        # type params
        self.turgor=0.01
        self.div_vol=100.0
        self.divaxes=0 # principal axis index
        
    def get_turgor_for_slider(self):
        return self.turgor*1000.0
    
    def set_turgor_from_slider(self,value):
        self.turgor=value/1000.0
                              
                              
                              
                              
                              
