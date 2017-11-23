from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.arrays import vbo
import math
import numpy
import numpy.linalg as la
import random

class GLGridRenderer:
    def __init__(self, sig, integ, rng=None, alpha=0.5):
        self.sig = sig
        self.integ = integ
        self.rng = rng
        self.alpha = alpha

        self.size = self.sig.gridSize
        self.dim = self.sig.gridDim[1:]
        self.len = self.dim[0]*self.dim[1]
        self.orig = [self.sig.gridOrig[i] - self.sig.gridSize[i]*0.5 for i in range(3)]
        self.texture = glGenTextures(1)
        self.imageData = numpy.zeros(self.dim)

        # Create byte array for texture data: OpenGL requires 2^n size and 
        # square?
        self.texDim = max([2**math.ceil(math.log(d,2)) for d in self.dim[0:2]])
        # The texture coordinates to fit the actual data
        self.tyMax = self.dim[0]/self.texDim
        self.txMax = self.dim[1]/self.texDim
        # The texture data, 3 bytes per (x,y,z)
        self.byteImageData = numpy.zeros((self.texDim,self.texDim)+(3,)).astype(numpy.uint8)

    def init_gl(self):
       pass
    
    def renderNames_gl(self):
       pass
        
    def render_gl(self, selection=None):
        signalLevel = self.integ.signalLevel.reshape(self.sig.gridDim)
        self.imageData = signalLevel[:,:,:,math.floor(self.dim[2]*0.5)] #numpy.log(self.sig.level[:,:,:,0]+1e-8)
        if self.rng:
            mx = self.rng[1]
            mn = self.rng[0]
        else:
            mx = numpy.max(self.imageData)
            mn = numpy.min(self.imageData)
        if mx>mn:
            scale = 255.0/(mx-mn)
        else:
            scale = 1.0
        #print "mx = %f, mn = %f, scale = %f"%(mx,mn,scale)
        self.imageData = numpy.clip((self.imageData - mn)*scale, 0.0, 255.0)
        #print "Signal grid range = %f to %f"%(mn,mx) 
        for s in range(self.sig.nSignals):
            self.byteImageData[0:self.dim[0],0:self.dim[1],s] = self.imageData[s,:,:].astype(numpy.uint8)

        glDisable(GL_CULL_FACE)
        glEnable(GL_TEXTURE_2D)
        glDisable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)
        glBindTexture(GL_TEXTURE_2D, self.texture)
        #Load the Data into the Texture
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, self.texDim, self.texDim, 0, GL_RGB, GL_UNSIGNED_BYTE, self.byteImageData )
#        glTexSubImage2Dub( GL_TEXTURE_2D, 0, 0,0, GL_RED, GL_UNSIGNED_BYTE, self.imageData )
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        #Define some Parameters for this Texture
        
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA ,GL_ONE_MINUS_SRC_ALPHA)
        glColor4f(1.0,1.0,1.0,self.alpha)
        glBegin(GL_QUADS)
        glTexCoord2f(0.0, 0.0)
        glVertex3f(self.orig[0], self.orig[1],  0.0)    # Bottom Left Of The Texture and Quad
        glTexCoord2f(0.0, self.tyMax)
        glVertex3f(self.orig[0]+self.dim[0]*self.size[0],  self.orig[1],  0.0)# Bottom Right Of The Texture and Quad
        glTexCoord2f(self.txMax, self.tyMax)
        glVertex3f( self.orig[0]+self.dim[0]*self.size[0], self.orig[1]+self.dim[1]*self.size[1],  0.0) # Top Right Of The Texture and Quad
        glTexCoord2f(self.txMax, 0.0)
        glVertex3f(self.orig[0] , self.orig[1]+self.dim[1]*self.size[1],  0.0) 
        glEnd()
        glEnable(GL_LIGHTING)
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_BLEND)
        glDisable(GL_TEXTURE_2D)

class GLPlantSignalRenderer:
    # Render signals as chanels rgb
    def __init__(self, sim, chanIdx):
        self.sim = sim
        self.wallcol = [0, 0, 0]
        self.nodecol = [0, 0, 0]
        self.cellcol = [1, 1, 1] 
	self.chanIdx = chanIdx
            
    def init_gl(self):
        pass    

    def renderNames_gl(self):
	pass

    def render_gl(self, selection=None):
	# Render signals[chanIdx] as rgb
	#
	maxSig = [0,0,0]
        for i in range(len(self.chanIdx)):
            for (cid,cellState) in self.sim.cellStates.items():
		s = cellState.signals[self.chanIdx[i]]
		if s > maxSig[i]:
		    maxSig[i] = s
	    if maxSig[i]==0:
		maxSig[i]=1.0

	for (cid,cellState) in self.sim.cellStates.items():
	    col = [0,0,0]
	    for i in range(len(self.chanIdx)):
		col[i] = cellState.signals[self.chanIdx[i]]/maxSig[i]
            glColor3fv(col)
            glBegin(GL_POLYGON)
            for p in cellState.nodep:
                glVertex2d(p[0], p[1])
            glEnd()
	    
class GLPlantRenderer:
    def __init__(self, sim):
        self.sim = sim
        self.wallcol = [0, 0, 0]
        self.nodecol = [0, 0, 0]
        self.cellcol = [1, 1, 1] 
            
    def init_gl(self):
        pass    

    def renderNames_gl(self):
	# Render cell_ids for selection
	#
        for (cid,cellState) in self.sim.cellStates.items():
            glColor3fv(self.cellcol)
            glPushName(cid)
            glBegin(GL_POLYGON)
            for p in cellState.nodep:
                glVertex2d(p[0], p[1])
            glEnd()
	    glPopName()

    def render_gl(self, selection=None):
        # Render model using PyOpenGL
	#
            # principle axes
            #scale=10.0
            #c=cell.get_centre()
            #glLineWidth(1)
            #glColor3f(1.0,0.0,0.0) #red
            #glBegin(GL_LINES)
            #glVertex2d(c[0],c[1])
            #glVertex2d(c[0]+cell.pa1[0]*scale,c[1]+cell.pa1[1]*scale)
            #glColor3f(0.0,0.0,1.0) #blue
            #glVertex2d(c[0],c[1])
            #glVertex2d(c[0]+cell.pa2[0]*scale,c[1]+cell.pa2[1]*scale)
            #glEnd() 
	glDisable(GL_LIGHTING)
	glDisable(GL_DEPTH_TEST)
	cellStates = self.sim.cellStates

	#for (cid,cellState) in cellStates.items():
	#    cellcol = cellState.color        
	#    glColor3fv(cellcol)
        #    glBegin(GL_POLYGON)
        #    for p in cellState.nodep:
        #        glVertex2d(p[0], p[1])
        #    glEnd()
	
	for (cid,cellState) in cellStates.items():
 	    col = self.wallcol
	    for (wp1,wp2) in cellState.wallp:
	        glColor3fv(col)
	        glLineWidth(2)
		glBegin(GL_LINES)
		glVertex2d(wp1[0], wp1[1])
		glVertex2d(wp2[0], wp2[1])
		glEnd()

	if selection>0 and cellStates.has_key(selection):
	    scell = cellStates[selection]
            glColor3fv([1,0,0])
	    for (wp1,wp2) in scell.wallp:
	        glLineWidth(2)
	        glBegin(GL_LINES)
	        glVertex2d(wp1[0], wp1[1])
	        glVertex2d(wp2[0], wp2[1])
	        glEnd()

	# Draw nodes as points
	#glColor3fv(self.nodecol)
        #glPointSize(3)
        
        #for node in self.model.nodes:
        #    p = node.get_pos()
        #    glBegin(GL_POINTS)
        #    glVertex2d(p[0], p[1])
        #    glEnd()
        #    r = node.get_rot()
        #    glBegin(GL_LINES)
        #    glVertex2d(p[0], p[1])
        #    glVertex2d(p[0] + math.cos(r), p[1] + math.sin(r))
        #    glEnd()
       
        glEnable(GL_LIGHTING)
       

    
    
class GLBacteriumRenderer:
    def __init__(self, sim, properties=None, scales = None):
        self.ncells_list = 0
        self.ncells_names_list = 0
        self.names_list_step = -1
        self.list_step = -1
        self.dlist = None
        self.dlist_names = None
        self.cellcol = [1, 1, 1] 
        self.sim = sim 
        self.quad = gluNewQuadric()
        self.properties = properties
        self.scales = scales

    def init_gl(self):
        pass

    def build_list(self, cells):
        if self.dlist:
            glDeleteLists(self.dlist, 1)
        index = glGenLists(1)
        glNewList(index, GL_COMPILE)
        self.render_cells()
        glEndList()
        self.dlist = index

    def build_list_names(self, cells):
        if self.dlist_names:
            glDeleteLists(self.dlist_names, 1)
        index = glGenLists(1)
        glNewList(index, GL_COMPILE)
        self.render_cell_names()
        glEndList()
        self.dlist_names = index

    def render_gl(self, selection=None):
        # Uncomment this to directly draw all cells each frame
        #self.render_cells()

        cells = self.sim.cellStates.values()
        if self.sim.stepNum > self.list_step:
            self.build_list(cells)
            #self.ncells_list = len(cells)
            self.list_step = self.sim.stepNum

        glCallList(self.dlist)
        
        if self.sim.cellStates.has_key(selection):
            self.render_selected_cell(selection)
        self.render_constraints()
	
    def render_constraints(self):
        phys = self.sim.phys
        glDisable(GL_LIGHTING)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LINE_SMOOTH)
        glLineWidth(1.0)
        #glBegin(GL_LINES)
        glColor3f(1,0,0)
        for i in range(phys.n_planes):
            p = phys.plane_pts[i]
            nn = phys.plane_norms[i]
            n = [nn[i] for i in range(3)]

            # Make an orthonormal basis in the plane
            v1 = numpy.cross([0,0,1],n)
            if la.norm(v1)<1e-6:
                v1 = numpy.cross([0,1,0],n)
            v2 = numpy.cross(v1,n)

            glBegin(GL_LINE_STRIP)
            glVertex3f(p[0], p[1], p[2])
            #glVertex3f(p[0]+n[0]*5, p[1]+n[1]*5, p[2]+n[2]*5)
            glVertex3f(p[0]+v1[0]*5, p[1]+v1[1]*5, p[2]+v1[2]*5)
            glVertex3f(p[0]+(v1[0]+v2[0])*5, p[1]+(v1[1]+v2[1])*5, p[2]+(v1[2]+v2[2])*5)
            glVertex3f(p[0]+v2[0]*5, p[1]+v2[1]*5, p[2]+v2[2]*5)
            glVertex3f(p[0], p[1], p[2])
            glEnd()
        #glEnd()
        glEnable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)


    def render_selected_cell(self, selection):
        glEnable(GL_DEPTH_TEST)
        cellcol = [1,0,0]
        cell = self.sim.cellStates[selection]
        l = cell.length
        r = cell.radius

        (e1,e2) = cell.ends
        ae1 = numpy.array(e1)
        ae2 = numpy.array(e2)
        zaxis = numpy.array([0,0,1])
        caxis = numpy.array(cell.dir) #(ae2-ae1)/l
        rotaxis = numpy.cross(caxis, zaxis)
        rotangle = numpy.arccos(numpy.dot(caxis,zaxis))
           
        # draw the outlines antialiased in RED
        glColor3f(1.0, 0.0, 0.0)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA ,GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_LINE_SMOOTH)
        glLineWidth(8.0)
        # draw wireframe for back facing polygons and cull front-facing ones
        glPolygonMode(GL_BACK, GL_FILL)
        glEnable(GL_CULL_FACE)
        glCullFace(GL_FRONT)
        glDepthFunc(GL_LEQUAL)

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glTranslatef(e1[0],e1[1],e1[2])
        gluSphere(self.quad, r, 8, 8)
        glRotatef(-rotangle*180.0/numpy.pi, rotaxis[0], rotaxis[1], rotaxis[2])
        gluCylinder(self.quad, r, r , l, 8, 1)
        glPopMatrix() 
        glPushMatrix()
        glTranslatef(e2[0],e2[1],e2[2])
        gluSphere(self.quad, r, 8, 8)
        glPopMatrix() 

        glDisable(GL_DEPTH_TEST)



    def renderNames_gl(self, selection=None):
        #self.render_cell_names()
        cells = self.sim.cellStates.values()
        if self.sim.stepNum > self.names_list_step:
            self.build_list_names(cells)
            self.names_list_step = self.sim.stepNum
        glCallList(self.dlist_names)

    def render_cell_names(self):
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
        glDepthFunc(GL_LESS)
        glDisable(GL_CULL_FACE)
        glPolygonMode(GL_FRONT, GL_FILL)
        glDisable(GL_LINE_SMOOTH)
        glDisable(GL_BLEND)

        for cell in self.sim.cellStates.values():
            l = cell.length
            r = cell.radius

            (e1,e2) = cell.ends
            ae1 = numpy.array(e1)
            ae2 = numpy.array(e2)
            zaxis = numpy.array([0,0,1])
            caxis = numpy.array(cell.dir) #(ae2-ae1)/l
            rotaxis = numpy.cross(caxis, zaxis)
            rotangle = numpy.arccos(numpy.dot(caxis,zaxis))
	    
            cid = cell.id
            glPushName(cid) 
	
            glMatrixMode(GL_MODELVIEW)
            glPushMatrix()
            glTranslatef(e1[0],e1[1],e1[2])
            glRotatef(-rotangle*180.0/numpy.pi, rotaxis[0], rotaxis[1], rotaxis[2])
            gluCylinder(self.quad, r, r , l, 8, 1)
            gluSphere(self.quad, r, 8, 8)
            glPopMatrix() 
            glPushMatrix()
            glTranslatef(e2[0],e2[1],e2[2])
            gluSphere(self.quad, r, 8, 8)
            glPopMatrix() 

	    glPopName()

        glEnable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)


    def render_cells(self, selection=None):
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)

        cells = self.sim.cellStates.values()
        for cell in cells:
            l = cell.length
            #r = cell.radius*2.0
            r = cell.radius

            (e1,e2) = cell.ends
            ae1 = numpy.array(e1)
            ae2 = numpy.array(e2)
            zaxis = numpy.array([0,0,1])
            caxis = numpy.array(cell.dir) #(ae2-ae1)/l
            rotaxis = numpy.cross(caxis, zaxis)
            rotangle = numpy.arccos(numpy.dot(caxis,zaxis))
           
            cid = cell.id
            cellcol = cell.color #self.cellcol #[random.uniform(0,1), random.uniform(0,1), random.uniform(0,1)] 
            if self.properties:
                cellcol = []
                for p in self.properties:
                    if hasattr(cell,p):
                        cellcol.append(getattr(cell,p))
                    else:
                        cellcol.append(0)
                for i in range(3):
                    cellcol[i] *= self.scales[i]
                    cellcol[i] = min(1,cellcol[i])
        

            # draw the outlines antialiased in black
            glColor3f(0.0, 0.0, 0.0)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA ,GL_ONE_MINUS_SRC_ALPHA)
            glEnable(GL_LINE_SMOOTH)
            glLineWidth(8.0)
            # draw wireframe for back facing polygons and cull front-facing ones
            glPolygonMode(GL_BACK, GL_FILL)
            glEnable(GL_CULL_FACE)
            glCullFace(GL_FRONT)
            glDepthFunc(GL_LEQUAL)

            glMatrixMode(GL_MODELVIEW)
            glPushMatrix()
            glTranslatef(e1[0],e1[1],e1[2])
            gluSphere(self.quad, r, 8, 8)
            glRotatef(-rotangle*180.0/numpy.pi, rotaxis[0], rotaxis[1], rotaxis[2])
            gluCylinder(self.quad, r, r , l, 8, 1)
            glPopMatrix() 
            glPushMatrix()
            glTranslatef(e2[0],e2[1],e2[2])
            gluSphere(self.quad, r, 8, 8)
            glPopMatrix() 

            glDepthFunc(GL_LESS)
            glDisable(GL_CULL_FACE)
            glPolygonMode(GL_FRONT, GL_FILL)
            glDisable(GL_LINE_SMOOTH)
            glDisable(GL_BLEND)

            glColor3fv(cellcol)
            glMatrixMode(GL_MODELVIEW)
            glPushMatrix()
            glTranslatef(e1[0],e1[1],e1[2])
            glScalef(0.8,0.8,0.8)
            gluSphere(self.quad, r, 8, 8)
            #glScalef(1.25,1.0,1.0)
            glRotatef(-rotangle*180.0/numpy.pi, rotaxis[0], rotaxis[1], rotaxis[2])
            gluCylinder(self.quad, r, r , l*1.25, 8, 1)
            glPopMatrix() 
            glPushMatrix()
            glTranslatef(e2[0],e2[1],e2[2])
            glScalef(0.8,0.8,0.8)
            gluSphere(self.quad, r, 8, 8)
            glPopMatrix() 
	
	# draw contact points
	if False: #hasattr(cells[0], 'contacts'):
	    glDisable(GL_DEPTH_TEST)
	    glDisable(GL_LIGHTING)
            for cell in cells:
		contacts = cell.contacts
		glBegin(GL_LINES)
		for ct in contacts:
		    glColor3fv(ct[6:9])
		    glVertex3fv(ct[0:3])
		    glVertex3fv(ct[3:6])
		glEnd()
		glBegin(GL_POINTS)
		for ct in contacts:
		    glColor3fv(ct[6:9])
		    glVertex3fv(ct[0:3])
		glEnd()
	    glEnable(GL_DEPTH_TEST)
	    glEnable(GL_LIGHTING)

            
class GLStaticMeshRenderer:
    def __init__(self, mesh, regul):
        self.mesh = mesh
        self.regul = regul
        self.vbo = None
        self.vboOffsets = {}
        self.nverts = 0

    def init_gl(self):
        self.makeVBOs()

    def makeVBOs(self):
        #print "Making VBO for triangle mesh..."
        cells = self.mesh.getCells()
        self.nverts = 0
        for c in cells:
            self.nverts += len(c.tris)*3
    
        self.normals = numpy.zeros((self.nverts,3))
        #ni = 0
        #for c in cells:
         #   for t in c.tris:
        #    v1,v2,v3 = self.mesh.verts[t[0]], self.mesh.verts[t[1]], self.mesh.verts[t[2]]
        #    v1 = numpy.array(v1)
        #    v2 = numpy.array(v2)
        #    v3 = numpy.array(v3)
        #    n = -numpy.cross((v2-v1),(v3-v1))
        #    n /= math.sqrt(numpy.dot(n,n))
        #    self.normals[ni,:] = n
        #    self.normals[ni+1,:] = n
        #    self.normals[ni+2,:] = n
        #    ni += 3
    
        varr = []
        voff = 0
        for c in cells:
            self.vboOffsets[c.getId()] = voff
            for t in c.tris:
                for ti in t:
                    for i in xrange(3):
                        varr.append(self.mesh.verts[ti][i])
                    for i in xrange(3):
                        varr.append(self.normals[voff][i])
                    for i in xrange(2):
                        varr.append(0.0) # pad to align to 32byte blocks
                    voff += 1
        self.vbo = vbo.VBO(numpy.array(varr,'f'))

    def render_gl(self, selection=None):
        if not self.vbo:
            self.init_gl()
        
        states = self.regul.cellStates
    
        glDisable(GL_LIGHTING)
        glEnable(GL_DEPTH_TEST)
        glShadeModel(GL_FLAT)
        glCullFace(GL_FRONT)
        glEnable(GL_CULL_FACE)
    
        self.vbo.bind()
        glEnableClientState(GL_VERTEX_ARRAY);
        #glEnableClientState(GL_NORMAL_ARRAY);
    
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
    
        glVertexPointer(3, GL_FLOAT, 32, self.vbo )
        #glNormalPointer(GL_FLOAT, 32, self.vbo+12 )
        for c in self.mesh.getCells():
            cid = c.getId()
            if selection!=cid:
                s =  states[cid]
		if len(s.signals)>0:
		    glColor4f(s.signals[0],0,0,1)
                else:
                    glColor4fv(c.color) #c.color[0:3]) 
            else:
                glColor4f(1,0,0,1)
            glDrawArrays(GL_TRIANGLES, self.vboOffsets[cid], len(c.tris)*3)
        
        glDisable(GL_POLYGON_OFFSET_FILL);
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glColor4f(0,0,0,0.5); 
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_TRIANGLES, 0, self.nverts)
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
        self.vbo.unbind()
        glDisableClientState(GL_VERTEX_ARRAY);
        #glDisableClientState(GL_NORMAL_ARRAY);
    
    #    glEnable(GL_LIGHTING)
            #glDisable(GL_BLEND) #DEPTH_TEST)
    #    glDisable(GL_DEPTH_TEST)
    
    def renderNames_gl(self):
        if not self.vbo:
            self.init_gl()
           
        glDisable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)
        #glCullFace(GL_FRONT)
        glDisable(GL_CULL_FACE)
        
        self.vbo.bind()
        glEnableClientState(GL_VERTEX_ARRAY);
        
        glVertexPointer(3, GL_FLOAT, 32, self.vbo )
        for c in self.mesh.getCells():
            cid = c.getId()
            glPushName(cid)
            glDrawArrays(GL_TRIANGLES, self.vboOffsets[cid], len(c.tris)*3)
            glPopName()
          
        self.vbo.unbind()
        glDisableClientState(GL_VERTEX_ARRAY);
        glEnable(GL_LIGHTING)




class GLCelBacteriumRenderer:
    def __init__(self, sim, properties=None, scales=None):
        self.cellcolor = [0.4, 0.6, 0.5]
        self.sim = sim
        self.properties = properties
        self.scales = scales
        self.init_gl()

    def init_gl(self):
        # display lists for drawing cells
        self.cylinder = self.build_cylinder(16)
        self.hemisphere = self.build_hemisphere(16, 16)
        # set up the 1d texture used for cel-shading 'lighting'
        self.texture = glGenTextures(1)
        glBindTexture(GL_TEXTURE_1D, self.texture)
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
        self.tex_img = [0.2]*(7*3) + [0.6]*(10*3) + [0.9]*(15*3)
        glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 32, 0, GL_RGB, GL_FLOAT, self.tex_img)
        glEnable(GL_TEXTURE_1D)
        glBindTexture(GL_TEXTURE_1D, self.texture)
        # we'll do our own 'lighting'
        glDisable(GL_LIGHTING)

    def build_cylinder(self, n_phis):
        """Return the id of a display list that draws a cylinder.

        The cylinder is along the x axis from -1/2 to 1/2.  It has a
        circular cross-section with radius 1.

        n_phis -- The number of rectangular panels with which to
        approximate the cylinder.
        """
        index = glGenLists(1)
        phis = [float(i)*2.0*numpy.pi/float(n_phis) for i in range(n_phis+1)]
        phi_pairs = zip(phis, phis[1:])
        glNewList(index, GL_COMPILE)
        glBegin(GL_QUADS)
        for phi1,phi2 in phi_pairs:
            dot1 = min(max(numpy.cos(phi1), 0.0), 1.0)
            dot2 = min(max(numpy.cos(phi2), 0.0), 1.0)
            glTexCoord1f(dot1)
            glVertex3f(-0.5, numpy.sin(phi1), numpy.cos(phi1))
            glTexCoord1f(dot1)
            glVertex3f(0.5, numpy.sin(phi1), numpy.cos(phi1))
            glTexCoord1f(dot2)
            glVertex3f(0.5, numpy.sin(phi2), numpy.cos(phi2))
            glTexCoord1f(dot2)
            glVertex3f(-0.5, numpy.sin(phi2), numpy.cos(phi2))
        glEnd()
        glEndList()
        return index

    def build_hemisphere(self, n_phis, n_thetas):
        """Return the id of a display list that draw a hemisphere.

        The hemisphere is in the +x half of the x axis and opens on the x=0 plane.  It has a radius of 1.

        n_phis -- Number of latitudes to use.  There are actually n_phis/2+1 latitudes.
        n_thetas -- Number of longitudes to use.
        """
        index = glGenLists(1)
        phis = [float(i)*numpy.pi/float(n_phis/2) for i in range(n_phis/2+1)]
        phi_pairs = zip(phis, phis[1:])
        thetas = [float(i)*numpy.pi/float(n_thetas) for i in range(-n_thetas/2, n_thetas/2+1)]
        theta_pairs = zip(thetas, thetas[1:])
        glNewList(index, GL_COMPILE)
        glBegin(GL_QUADS)
        for phi1,phi2 in phi_pairs:
            dot1 = min(max(numpy.cos(phi1), 0.0), 1.0)
            dot2 = min(max(numpy.cos(phi2), 0.0), 1.0)
            for th1,th2 in theta_pairs:
                glTexCoord1f(dot1)
                glVertex3f(numpy.sin(phi1)*numpy.cos(th1), numpy.sin(phi1)*numpy.sin(th1), numpy.cos(phi1))
                glTexCoord1f(dot2)
                glVertex3f(numpy.sin(phi2)*numpy.cos(th1), numpy.sin(phi2)*numpy.sin(th1), numpy.cos(phi2))
                glTexCoord1f(dot2)
                glVertex3f(numpy.sin(phi2)*numpy.cos(th2), numpy.sin(phi2)*numpy.sin(th2), numpy.cos(phi2))
                glTexCoord1f(dot1)
                glVertex3f(numpy.sin(phi1)*numpy.cos(th2), numpy.sin(phi1)*numpy.sin(th2), numpy.cos(phi1))
        glEnd()
        glEndList()
        return index

    def render_capsule(self, l, r):
        """Draw a capsule at the origin along the x axis.

        l -- length of the capsule.
        r -- radius of the capsule.
        """
        # draw cylinder
        glPushMatrix()
        glScalef(l, r, r)
        glCallList(self.cylinder)
        glPopMatrix()
        # draw +x hemisphere
        glPushMatrix()
        glTranslatef(l/2.0, 0, 0)
        glScalef(r, r, r)
        glCallList(self.hemisphere)
        glPopMatrix()
        # draw -x hemisphere
        glPushMatrix()
        glRotatef(180.0, 0, 0, 1)
        glTranslatef(l/2.0, 0, 0)
        glScalef(r, r, r)
        glCallList(self.hemisphere)
        glPopMatrix()

    def render_cell(self, cell, selection=None):
        l = cell.length
        r = cell.radius
        (e1,e2) = map(numpy.array, cell.ends)
        ori = (e1+e2)/2.0 # origin in the cell's coordinates
        z = numpy.array([0, 0, 1])
        axis = (e2-e1)/l # x axis in the cell's coordinates
        rotangle = numpy.arctan2(axis[1], axis[0]) # angle to rotate around the z axis

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()

        # move to the cell's coordinates
        glTranslatef(ori[0], ori[1], ori[2])
        glRotatef(rotangle*180.0/numpy.pi, z[0], z[1], z[2])

        # draw the cell in color
        color = getattr(cell, 'color', self.cellcolor)
        glColor3fv(color)
        self.render_capsule(l, r)

        # draw the outlines antialiased in black
        glColor3f(0.0, 0.0, 0.0)
        glLineWidth(2.0)
        if cell.id == selection:
            glColor3f(1.0, 0.3, 0.3)
            glLineWidth(3.0)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA ,GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_LINE_SMOOTH)
        # draw wireframe for back facing polygons and cull front-facing ones
        glPolygonMode(GL_BACK, GL_LINE)
        glEnable(GL_CULL_FACE)
        glCullFace(GL_FRONT)
        glDepthFunc(GL_LEQUAL)
        self.render_capsule(l, r)
        glDepthFunc(GL_LESS)
        glDisable(GL_CULL_FACE)
        glPolygonMode(GL_FRONT, GL_FILL)
        glDisable(GL_LINE_SMOOTH)
        glDisable(GL_BLEND)

        glPopMatrix()

    def render_cell_name(self, cell):
        l = cell.length
        r = cell.radius
        (e1,e2) = map(numpy.array, cell.ends)
        ori = (e1+e2)/2.0 # origin in the cell's coordinates
        z = numpy.array([0, 0, 1])
        axis = (e2-e1)/l # x axis in the cell's coordinates
        rotangle = numpy.arctan2(axis[1], axis[0]) # angle to rotate around the z axis

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()

        # move to the cell's coordinates
        glTranslatef(ori[0], ori[1], ori[2])
        glRotatef(rotangle*180.0/numpy.pi, z[0], z[1], z[2])

        # draw the cell in color
        glPushName(cell.id)
        self.render_capsule(l, r)
        glPopName()

        glPopMatrix()


    def render_gl(self, selection=None):
  	cells = self.sim.cellStates.values()
        for cell in cells: self.render_cell(cell, selection)

    def renderNames_gl(self, selection=None):
        cells = self.sim.cellStates.values()
        for cell in cells: self.render_cell_name(cell)







class GL2DBacteriumRenderer:

    def __init__(self, sim, properties=None, scales=None):
        self.ncells_list = 0
        self.ncells_names_list = 0
        self.dlist = None
        self.dlist_names = None
        self.sim = sim
        self.properties = properties
        self.scales = scales
        self.cellcolor = [0.4, 0.6, 0.5]
        self.circ_pts = [(math.cos(math.radians(th)), math.sin(math.radians(th))) for th in range(-80,90,20)]
        self.semi = self.build_semicircle()


    def build_semicircle(self):
        index = glGenLists(1)
        glNewList(index, GL_COMPILE)
        glVertex3f(0.5, -1, 0)
        for x,y in self.circ_pts:
            glVertex3f(0.5 + x, y, 0.0)
        glVertex3f(0.5, 1, 0)
        glEndList()
        return index

    def draw_cell(self, p, d, l, r):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        ang = math.atan2(d[1], d[0]) * 360.0 / (2.0*3.141593)
        glTranslatef(p[0], p[1], 0.0)
        glRotatef(ang, 0.0, 0.0, 1.0)
        glBegin(GL_POLYGON)
        glVertex3f(-l/2.0, -r, 0)
        glVertex3f(l/2.0, -r, 0)
        for x,y in self.circ_pts:
            glVertex3f(l/2.0 + x*r, y*r, 0.0)
        glVertex3f(l/2.0, r, 0)
        glVertex3f(-l/2.0, r, 0)
        for x,y in self.circ_pts:
            glVertex3f(-l/2.0 -x*r, -y*r, 0.0)
        glEnd()
        glPopMatrix()

    def build_list(self, cells):
        if self.dlist:
            glDeleteLists(self.dlist, 1)
        index = glGenLists(1)
        glNewList(index, GL_COMPILE)
        for cell in cells: self.render_cell(cell)
        glEndList()
        self.dlist = index

    def build_list_names(self, cells):
        if self.dlist_names:
            glDeleteLists(self.dlist_names, 1)
        index = glGenLists(1)
        glNewList(index, GL_COMPILE)
        for cell in cells: self.render_cell_name(cell)
        glEndList()
        self.dlist_names = index

    def render_cell(self, cell, selection=None):
        p,d,l,r = (cell.pos, cell.dir, cell.length, cell.radius)
        glMatrixMode(GL_MODELVIEW)
        glEnable(GL_DEPTH_TEST)
        glPushMatrix()
        glPolygonMode(GL_FRONT, GL_FILL)
        glColor3f(0.0, 0.0, 0.0)
        if selection == cell.id:
            glColor3f(0.2, 0.1, 0.8)
        glTranslatef(0.0, 0.0, 0.025)
        self.draw_cell(p, d, l, r)
        color = getattr(cell, 'color', [0.4, 0.4, 0.4])
        glColor3fv(color)
        glTranslatef(0.0, 0.0, 0.025)
        self.draw_cell(p, d, l, r-0.15)
        glPopMatrix()
        glDisable(GL_DEPTH_TEST)


    def render_cell_name(self, cell, selection=None):
        p,d,l,r = (cell.pos, cell.dir, cell.length, cell.radius)
        glPushName(cell.id)
        glEnable(GL_DEPTH_TEST)
        self.draw_cell(p, d, l, r)
        glPopName()
        glDisable(GL_DEPTH_TEST)


    def render_gl(self, selection=None):
        glDisable(GL_LIGHTING)
        glEnable(GL_POLYGON_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        cells = self.sim.cellStates.values()
        if len(cells)!=self.ncells_list:
            self.build_list(cells)
            self.ncells_list = len(cells)
        glCallList(self.dlist)
        #for cell in cells: self.render_cell(cell, selection)


    def renderNames_gl(self, selection=None):
        cells = self.sim.cellStates.values()
        if len(cells)!=self.ncells_names_list:
            self.build_list_names(cells)
            self.ncells_names_list = len(cells)
        glCallList(self.dlist_names)
        #for cell in cells: self.render_cell_name(cell, selection)



