import sys
import os
import math
import numpy as np
import pickle
import CellModeller
import subprocess
import string
import shutil
from CellModeller.Simulator import Simulator
import networkx
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib import units
from reportlab.lib.colors import Color
from reportlab.graphics.shapes import Circle

class CellModellerPDFGenerator(Canvas):
    # ---
    # Class that extends reportlab pdf canvas to draw CellModeller simulations
    # ---
    def __init__(self, name, network, bg_color):
        self.name = name
        self.data = data
        self.bg_color = bg_color
        Canvas.__init__(self, name)
        self.network = network
    
    def setup_canvas(self, name, world, page, page_center):
        worldx,worldy = world
        pagex,pagey = page
        page_centx,page_centy = page_center
        self.setPageSize((pagex*units.cm, pagey*units.cm))
        self.setFillColor(self.bg_color)
        self.rect(0, 0, pagex*units.cm, pagey*units.cm, fill=1)
        self.translate(pagex*units.cm/2.0, pagey*units.cm/2.0)
        self.translate(page_centx*units.cm, page_centy*units.cm)
        self.scale(float(pagex*units.cm)/worldx, float(pagey*units.cm)/worldy)

    def drawNode(self,pos,col):
        Canvas.setFillColor(self,col)
        Canvas.circle(self,pos[0],pos[1],0.3, fill = 1)

    def drawEdge(self,pos1,pos2):
        path = self.beginPath()
        path.moveTo(pos1[0],pos1[1])
        path.lineTo(pos2[0],pos2[1])
        self.drawPath(path, fill=1)

    def draw_frame(self, name, world, page, center):
        self.setup_canvas(name, world, page, center)
        self.draw_graph()
        self.showPage()
        self.save()

    def draw_graph(self):
        self.saveState()
        self.setLineWidth(0.003*units.cm)
        self.setStrokeColor((0.0,0.0,0.0))
        positions = list(networkx.get_node_attributes(G,'pos').values())
        colors = list(networkx.get_node_attributes(G,'color').values())
        for u,v,di in self.network.edges_iter(data=True):
            self.drawEdge(positions[u],positions[v])
        for n in self.network.nodes_iter():
            self.drawNode(positions[n],colors[n])
        self.restoreState()

def generate_network(G,num_cells,ct_tos,n_cts,pos,celldata):
    for i in range(0, num_cells): #make all nodes
        G.add_node(i, pos = pos[i], type=celldata[i], color = (1,0,0) if celldata[i]==0 else (0,1,0))
    for i in range(0, num_cells): #make all edges (when all nodes are present)
        for j in range(0, n_cts[i]):  #make all edges from contact_tos
            G.add_edge(i,ct_tos[i,j], width=8 , color = 'black')

def get_current_contacts(G, data):
    cs = data['cellStates']
    it = iter(cs)
    n = len(cs)
    cell_type={}
    pos_dict={}
    for it in cs:
        cell_type[cs[it].idx] = cs[it].cellType
        pos_dict[cs[it].idx] = cs[it].pos[0:2]
    modname = data['moduleName']
    moduleStr = data['moduleStr']
    sim = Simulator(modname, 0.0, moduleStr=moduleStr, saveOutput=False)
    sim.loadFromPickle(data)
    sim.phys.update_grid() # we assume local cell_centers is current
    sim.phys.bin_cells()
    sim.phys.cell_sqs = sim.phys.cell_sqs_dev.get() # get updated cell sqs
    sim.phys.sort_cells()
    sim.phys.sorted_ids_dev.set(sim.phys.sorted_ids) # push changes to the device
    sim.phys.sq_inds_dev.set(sim.phys.sq_inds)
    sim.phys.find_contacts(predict=False)
    sim.phys.get_cts()
    ct_pts=sim.phys.ct_pts #these are the points on the cell surface - they can be transformed into the global coordinate system
    ct_tos=sim.phys.ct_tos #this is the list of *some* contacts from each cell (only for lower cell_ids)
    ct_dists=sim.phys.ct_dists
    cell_cts=sim.phys.cell_n_cts #not really all the contacts of the cell, because the tos are only defined between a cell and ones with lower ids
    generate_network(G,n,ct_tos,cell_cts,pos_dict,cell_type)


bg_color = Color(1.0,1.0,1.0,alpha=1.0)


G = networkx.Graph()
fname = sys.argv[1]
data = pickle.load(open(fname,'rb'))
cs = data['cellStates']
it = iter(cs)
n = len(cs)
oname = fname.replace('.pickle','_graph.pdf')

print(("num_cells = "+str(n)))

cell_type={}
pos_dict={}
for it in cs:
    cell_type[cs[it].idx] = cs[it].cellType
    pos_dict[cs[it].idx] = cs[it].pos[0:2]

get_current_contacts(G, data)

print(("num_contacts = " + str(networkx.number_of_edges(G))))
degrees = list(G.degree().values())
print(("mean degree = " + str(np.mean(degrees))))

if list(networkx.get_edge_attributes(G,'color').items()):
    edges,ecolor = list(zip(*list(networkx.get_edge_attributes(G,'color').items())))
ncolor = list(networkx.get_node_attributes(G,'color').values())

pdf = CellModellerPDFGenerator(oname, G, bg_color)
world = (120,120)
page = (50,50)
center = (0,0)

pdf.draw_frame(oname, world, page, center)


#from this, I want the spatial position of every contact, and the graph of the cells that each cell is touching
#This can be all in one data structure - a graph where each vertex has a position, and each node has a cellState

