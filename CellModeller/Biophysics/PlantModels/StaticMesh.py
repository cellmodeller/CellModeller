import re
import numpy
import math
import random
import copy


class StaticCell:
	def __init__(self):
		self.id = -1
		self.nbs = []
		self.tris = []
		self.position = [0,0,0]
		self.vol = 0
		self.color = [random.uniform(0,1), random.uniform(0,1), random.uniform(0,1), 1]
		
	def pos(self):
		return self.position
	
	def volume(self):
		return self.vol
	
	def getId(self):
		return self.id
		
	def setGrowthRate(self, rate):
	    pass

class StaticMesh:

	def __init__(self, fname):
		self.verts = []
		self.tris = []
		self.cells = {}
		self.nbs = {}
		self.loadMesh(fname)
		del self.cells[-1]
		self.computeGeometry()
		for c in self.cells.values():
			self.nbs[c.id]= c.nbs

	def update(self, dt):
		pass
	
	def getCells(self):
		return self.cells.values()
	
	def getNeighbours(self):
		return self.nbs
	
	def divide(self, cell, asymm0, asymm1):
		pass

	def computeGeometry(self):
		print "Computing mesh geometry..."


		for c in self.cells.values():
			c.vol = 0
			c.area = 0
			c.position = numpy.zeros(3)
			for t in c.tris:
				v1,v2,v3 = self.verts[t[0]], self.verts[t[1]], self.verts[t[2]]
				v1 = numpy.array(v1)
				v2 = numpy.array(v2)
				v3 = numpy.array(v3)
				sA = numpy.cross((v2-v1),(v3-v2))
				A = math.sqrt(numpy.dot(sA,sA))*0.5
				c.vol += numpy.dot(sA,(v1+v2+v3))/18.0
				c.area += A
				c.position += (v1+v2+v3)/3.0
			c.position /= len(c.tris) #c.area 
			c.position = c.position.tolist()
			#print "Cell %d, position = %f,%f,%f, volume = %f, area = %f"%(c.position[0], c.position[1], c.position[2], c.getId(),c.vol,c.area)
				
	def loadMesh(self, file):		
		print "Loading Amira mesh file, %s" %(file)
		readVerts = False
		readTris = False
		readPatches = False

		for line in open(file):
			words = line.split()
		
			if readVerts:
				if len(self.verts)<numVerts:
					v = [float(s) for s in words[0:3]]
					self.verts.append(v)
			if readTris:
				if '}' in words or '{' in words:
					readTris = False
					innerCell = None
					outerCell = None
				else:
					if readPatches and innerCell and outerCell:
						#if innerCell==-1:
					#		innerCell = outerCell
					#	if outerCell==-1:
					#		outerCell = innerCell
						t = [int(s)-1 for s in words[0:3]]
						self.cells[innerCell].tris.append(t)
						tr = copy.copy(t)
						tr.reverse()
						self.cells[outerCell].tris.append(tr)
						self.tris.append(t)
			if readPatches:
				if 'InnerRegion' in words:
					n = re.findall('\d+',words[1])
					if n:
						innerCell = int(n[0])-1
					else:
						innerCell = -1
					if not self.cells.has_key(innerCell) and innerCell!=-1:
						self.cells[innerCell] = StaticCell()
						self.cells[innerCell].id = innerCell
				if 'OuterRegion' in words:
					n = re.findall('\d+',words[1])
					if n:
						outerCell = int(n[0])-1
					else:
						outerCell = -1
					if not self.cells.has_key(outerCell):
						self.cells[outerCell] = StaticCell()
						self.cells[outerCell].id = outerCell
					if innerCell and outerCell!=-1:
						self.cells[innerCell].nbs.append(outerCell)
					if outerCell and innerCell!=-1:
						self.cells[outerCell].nbs.append(innerCell)
					
		
			if 'Vertices' in words:
				readVerts = True
				readTris = False
				readPatches = False
				numVerts = int(words[1])
			elif 'Triangles' in words:
				readTris = True
				readVerts = False
				numTris = int(words[1])
			elif 'Patches' in words:
				readPatches = True
				readTris = False
				readVerts = False
				numPatches = int(words[1])
				innerCell=None
				outerCell=None
		
			
