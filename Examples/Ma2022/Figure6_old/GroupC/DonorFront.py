import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add

max_cells = 20000

delta = 2.0
delta_sig = 0.45

#Specify parameter for solving diffusion dynamics #Add
grid_dim = (120, 120, 12) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-240, -240, -24) # where to place the diffusion space onto simulation space

n_signals = 1
n_species = 1


def setup(sim):
    # Set biophysics module
    biophys = CLBacterium(sim, jitter_z=False, max_cells=max_cells, cgs_tol=1E-5,compNeighbours=True)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [10.0, 10.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)

    # Specify the initial cell and its location in the simulation
    dist = 5 #distance between cells on the grid
    radius = 60

    ################ Place producer at frontier: in total 437 cells; 219 recipients, 218 donors. almost 1:1 ratio.

    ###### To inoculate recipients behind frontier.
    coord = [] ### give an empty list for all possible coordinates.
    for m in range(-19,20):
        for n in range(-19,20):
            R = numpy.sqrt((m*dist)**2 + (n*dist)**2)
            if R < 55:
                points = (m,n)
                coord.append(points)
    points_rd = random.sample(coord,219) #### randomly pick 219 coordinates for recipients (producer)
    donors_within = list(set(coord) - set(points_rd)) ###### get complementary coordinates for donors behind the frontier.

    for element in points_rd:
        x = element[0]
        y = element[1]
        sim.addCell(cellType = 1, pos = (x*dist,y*dist,0) , dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4) ## inoculate recipients

    ### To inoculate donors at frontier.
    for i in range(-19,20):
        for j in range(-19,20):
            R = numpy.sqrt((i*dist)**2 + (j*dist)**2)
            if R >=55 and R < radius:
                sim.addCell(cellType=0, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4) ## inoculate donors at front

    for element_donor in donors_within:
        e = element_donor[0]
        f = element_donor[1]
        sim.addCell(cellType=0, pos=(e*dist,f*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4) ## inoculate donors behind.



            #elif R < 18:
            #    sim.addCell(cellType=1, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4)


     #producer at frontier: in total 56 donors; and 249 recipients, to make it roughly 1:1, need to add 97 more donors.



    #ACCEPTORS:
    #sim.addCell(cellType=0, pos=(-5,-5,0), dir=(1,0,0), length=delta,rad=0.4)
    #sim.addCell(cellType=0, pos=(5,5,0), dir=(1,0,0), length=delta,rad=0.4)
    #DONORS:
    #sim.addCell(cellType=1, pos=(5,-5,0), dir=(1,0,0), length=delta,rad=0.4)
    #sim.addCell(cellType=1, pos=(-5,5,0), dir=(1,0,0), length=delta,rad=0.4)

    # Add some objects to draw the models
    #therenderer = Renderers.GLBacteriumRenderer(sim)
    #sim.addRenderer(therenderer)
    #sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    #sim.addRenderer(sigrend) #Add

    # Specify how often data is saved
    sim.pickleSteps = 10

    sim.countHGT = True #Turns on HGT event counting
    sim.HGTevents = 0 #You need to start the counter here


def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell.length + random.gauss(delta,delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    #if cell.cellType == 0:
    #    cell.color = (1.0, 0.1, 0.1) #red acceptor, grows faster
    #    cell.growthRate = 1.1
    if cell.cellType == 0:
        cell.color = (0.98, 1.0, 0.15)#yellow consumer, donor
    elif cell.cellType == 1:
        cell.color = (0, 0.2, 0.7) #blue producer, recipient

    # Specify initial concentration of chemical species
    cell.species[:] = [0.0]*n_species
    # Specify initial concentration of signaling molecules
    cell.signals[:] = [0.0]*n_signals

cl_prefix = \
    '''
        const float Da = 1.0f;
        const float ka = 5.0f;
        float  alpha_in = species[0];
        float  alpha = signals[0];
        '''

# Da = diffusion rate of alpha through the cell membrane

def specRateCL(): # Chemical only produced in cellType 1
    global cl_prefix
    return cl_prefix + '''
        if (cellType==0){
        rates[0] = Da*(alpha-alpha_in)*area/gridVolume;
        } else {
        rates[0] = ka + Da*(alpha-alpha_in)*area/gridVolume;
        }
        '''

def sigRateCL(): #Add
    global cl_prefix
    return cl_prefix + '''
        rates[0] = -Da*(alpha-alpha_in)*area/gridVolume;
        '''

def update(cells,sim): #The module regulator now passes the simulator in ONLY IF countHGT is True
    v_max = 0.9
    Km = 0.1 #This can be thought of as the degree of mutualism - turns mutualism off when set to 0.

    #Iterate through each cell and flag cells that reach target size for division
    #for (id, cell) in cells.items():
    #cell.color = [1,1,1] #[0.1+cell.species[0]/3.0, 0.1+cell.species[1]/3.0, 0.1]
    #    if cell.cellType==0: #consumer growth rate depends on
    #        cell.growthRate = 0.5 + v_max * cell.species[0] / (Km + cell.species[0])
    #    if cell.volume > cell.targetVol:
    #        cell.divideFlag = True

    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
        infect_chances = 0
        if cell.cellType == 0:
            cell.growthRate = v_max * cell.species[0] / (Km + cell.species[0])
        if cell.cellType == 1: #if a cell is an acceptor producer,
            cell.growthRate = 1.0
            for index in cell.neighbours:#loop through all contacts
                if cells[index].cellType == 0: #if donor or transconjugant is in contact
                    if random.random() < 0.0005: # constant probability per unit time
                        sim.HGTevents += 1
                        cell.cellType = 2 #become transconjugant cyan
                        cell.growthRate = 1.0
                        cell.color = (0.32, 1.0, 0.96) ## cyan transconjugants

        if sim.stepNum > 400: #antibiotic added at time 400
            if cell.cellType == 1: #if a cell is an recipient,
                cell.growthRate = 0.0
            elif cell.cellType == 2: #if a cell is a transconjugant -- producer:
                cell.growthRate = 1
            else:
                cell.growthRate = v_max * cell.species[0] / (Km + cell.species[0])

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)
