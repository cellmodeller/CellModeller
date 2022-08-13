
def simulation_script(script_name, growth_rate):

    script = """
import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add


max_cells = 10000

#Specify parameter for solving diffusion dynamics #Add
grid_dim = (80, 80, 8) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-160, -160, -16) # where to place the diffusion space onto simulation space

n_signals = 1
n_species = 2


def setup(sim):
    sim.dt = 0.016
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, gamma=50, max_cells=max_cells, cgs_tol=1E-5,)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [60.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)
    
    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)
    
    # Specify the initial cell and its location in the simulation
    # RFP
    sim.addCell(cellType = 0,
                    pos = (50.5681509,-52.39245283, 0),
                    length = 3.988726,
                    dir = (-0.704020057,0.71018009,0),
                    rad = 0.547956)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (-49.86211765,-52.47388235, 0),
                    length = 2.042367535,
                    dir = (-0.964712925,0.263303954,0),
                    rad = 0.441133915)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (16.1188,-32.4084, 0),
                    length = 3.087266144,
                    dir = (-0.792682284,0.60963497,0),
                    rad = 0.521352221)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (18.60980645,-34.39703226, 0),
                    length = 3.261936006,
                    dir = (-0.881637656,0.471926947,0),
                    rad = 0.506656969)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (21.701589,-11.91791781, 0),
                    length = 3.696132851,
                    dir=(-0.368698908,0.929548877,0),
                    rad = 0.538698307)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (-52.89227184,-7.42741748, 0),
                    length = 2.918020639,
                    dir = (-0.518911374,0.854828045,0),
                    rad = 0.466708462)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (-50.01085217,-4.03575652, 0),
                    length = 2.825037476,
                    dir = (-0.254051842,0.967190602,0),
                    rad = 0.541402802)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (26.6267826,10.55304348, 0),
                    length = 1.819007216,
                    dir = (-0.955930987,0.293591464,0),
                    rad = 0.500360774)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (-12.43705882,12.448, 0),
                    length = 3.242584654,
                    dir = (-0.489771038,0.871851094,0),
                    rad = 0.559892112)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (15.71698701,29.6718961, 0),
                    length = 1.971501085,
                    dir = (-0.830593054,0.55687986,0),
                    rad = 0.517248043)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (11.22030769,38.9624615, 0),
                    length = 4.474477554,
                    dir = (-0.510736819,0.859737112,0),
                    rad = 0.472865773)  #Add RFP
    sim.addCell(cellType = 0,
                    pos = (-22.467,40.3225, 0),
                    length = 2.531593497,
                    dir = (-0.472729053,0.881207832,0),
                    rad = 0.508269211)  #Add RFP

        
    #YFP
    sim.addCell(cellType = 1,
                    pos = (55.276866,-80.360742268, 0),
                    length = 3.988726066,
                    dir= (-0.704020057,0.71018009,0),
                    rad = 0.547956075) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (-7.76806015,-79.58724812, 0),
                    length = 2.042367535,
                    dir= (-0.964712925,0.263303954,0),
                    rad = 0.441133915) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (-4.53821538,-77.036184615, 0),
                    length = 3.087266144,
                    dir= (-0.792682284,0.60963497,0),
                    rad = 0.521352221) #Add   #YFP                      
    sim.addCell(cellType = 1,
                    pos = (6.6477346,-53.33664455, 0),
                    length = 3.261936006,
                    dir= (-0.881637656,0.471926947,0),
                    rad = 0.506656969) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (10.58742857,-60.12647619, 0),
                    length = 3.696132851,
                    dir= (-0.368698908,0.929548877,0),
                    rad = 0.538698307) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (24.7379592,-41.85273469, 0),
                    length = 2.918020639,
                    dir= (-0.518911374,0.854828045,0),
                    rad = 0.466708462) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (25.4927059,-39.40254118, 0),
                    length = 2.825037476,
                    dir= (-0.254051842,0.967190602,0),
                    rad = 0.541402802) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (52.5072239,-34.82149254, 0),
                    length = 1.819007216,
                    dir= (-0.955930987,0.293591464,0),
                    rad = 0.500360774) #Add   #YFP
    sim.addCell(cellType = 1,
                    pos = (54.8251429,-30.03428571, 0),
                    length = 3.242584654,
                    dir= (-0.489771038,0.871851094,0),
                    rad = 0.559892112) #Add   #YFP                    
    sim.addCell(cellType = 1,
                    pos = (115.669,-13.68671428, 0),
                    length = 3.626993993,
                    dir= (0.015134651,0.999885465,0),
                    rad = 0.418235446) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (84.22063163,33.17894737, 0),
                    length = 2.037328247,
                    dir= (0.622267577,0.782804613,0),
                    rad = 0.371262294) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (90.95846147,35.97261538, 0),
                    length = 2.112031993,
                    dir= (-0.973457276,0.228868808,0),
                    rad = 0.33794242) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (85.67386206,37.24331035, 0),
                    length = 6.4043063,
                    dir= (-0.953934113,0.300016181,0),
                    rad = 0.511947448) #Add   #YFP  
    sim.addCell(cellType = 1,
                    pos = (34.41949579,40.74218488, 0),
                    length = 2.992585036,
                    dir= (0.350178972,0.936682811,0),
                    rad = 0.534955005) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (-14.24070588,43.23752942, 0),
                    length = 7.450685974,
                    dir= (0.726143354,0.687543329,0),
                    rad = 0.545179775) #Add   #YFP  
    sim.addCell(cellType = 1,
                    pos = (32.0521081,46.80821621, 0),
                    length = 5.835986471,
                    dir= (0.233742704,0.972298487,0),
                    rad = 0.519294116) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (-59.0077647,51.3957647, 0),
                    length = 2.765778018,
                    dir= (0.871639353,0.490147772,0),
                    rad = 0.325392278) #Add   #YFP   
    sim.addCell(cellType = 1,
                    pos = (-58.5728,74.68701544, 0),
                    length = 1.741159252,
                    dir= (0.132871888,0.991133221,0),
                    rad = 0.499394094) #Add   #YFP 
    sim.addCell(cellType = 1,
                    pos = (-58.56826966,76.51182016, 0),
                    length = 2.552424542,
                    dir= (-0.042043577,0.999115778,0),
                    rad = 0.478532423) #Add   #YFP     
    sim.addCell(cellType = 1,
                    pos = (105.516536,76.80115461, 0),
                    length = 2.56366508,
                    dir= (0.912424336,0.409245442,0),
                    rad = 0.501524163) #Add   #YFP 


    
    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    sim.addRenderer(sigrend) #Add
    
    sim.pickleSteps = 200
    sim.saveOutput=True

def init(cell):
    # Specify mean and distribution of initial cell size
    # 0: RFP(Attacker)  , 1: YFP (susceptible)
    if cell.cellType==0:
	    cell.targetVol = random.gauss(5.19,3.32)
	    # Specify growth rate of cells
	    cell.growthRate = random.gauss(0.03, 0.08)
    elif cell.cellType==1:
	    cell.targetVol = numpy.random.lognormal(1.47,0.62)
	    # Specify growth rate of cells
	    cell.growthRate = """ + str(growth_rate) + """
    #color
    if cell.cellType == 0: cell.color = (1.0,0.0,0.0) #RFP
    elif cell.cellType == 1: cell.color = (1.0,1.0,0.0) #YFP   
    # Specify initial concentration of chemical species
    cell.species[:] = [0.0]*n_species
    # Specify initial concentration of signaling molecules
    cell.signals[:] = [0.0]*n_signals




#Da: diffusion rates
# Ka,Kb: production rates
#alpha_in , alpha: chemical concentrations both outside and inside the cells.
cl_prefix = \
    '''
        const float Da = 1.f;
        const float ka = 1.f;
        const float kb = 1.f;
        
        float  alpha_in = species[0];
        float  alpha = signals[0];
        
        float beta_in = species[1];
        '''


# Da = diffusion rate of alpha through the cell membrane
# intracellular chemical reactions

def specRateCL(): # Add if/else, new species
    global cl_prefix
    return cl_prefix + '''
        rates[0] = ka + Da*(alpha_in)*area/gridVolume;
        rates[1] = kb*(beta_in)*(alpha-alpha_in)*area/gridVolume;

        '''


# behaviour of the extracellular signals.
def sigRateCL(): #Add
    global cl_prefix
    return cl_prefix + '''
    	if (cellType==0){
        rates[0] =  -Da*(alpha-alpha_in)*area/gridVolume;
        }else{
        rates[0]=0;
        }
        
        
        '''


def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    signalThreshold=0.0005
    for (id, cell) in cells.items():
        if cell.cellType==0:
            cell.growthRate = random.gauss(0.03 ,0.08)
            if cell.growthRate <=0.03:
                cell.species[0]=0
                cell.signals[0]=0
            if cell.volume > cell.targetVol:
                cell.divideFlag = True
                
        elif cell.cellType==1:
            cell.growthRate =  """ + str(growth_rate) + """
            if cell.volume > cell.targetVol:
                cell.divideFlag = True
            if cell.signals[0]>=signalThreshold:
                cell.cellType=2    
                
                           
        elif cell.cellType==2:
            cell.growthRate = 0

        if cell.cellType == 0: cell.color = (1.0,0.0,0.0) #RFP
        elif cell.cellType == 1: cell.color = (1.0,1.0,0.0) #YFP
        elif cell.cellType == 2: cell.color = (0.0,0.0,1.0) #Blue 

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    # 0: RFP(Attacker)  , 1: YFP (susceptible)
    if d1.cellType==0:
	    d1.targetVol = random.gauss(5.19,3.32)
	    d1.rad = random.gauss(0.57,0.31)
	    d1.growthRate = random.gauss(0.03 ,0.08)
    elif d1.cellType==1:
	    d1.targetVol = numpy.random.lognormal(1.47,0.62)
	    d1.rad = numpy.random.lognormal(-1.2,0.72)	
	    d1.growthRate = """ + str(growth_rate) + """
	    
    if d2.cellType==0:
	    d2.targetVol = random.gauss(5.19,3.32)
	    d2.rad = random.gauss(0.57,0.31)
	    d2.growthRate = random.gauss(0.03 ,0.08)
    elif d2.cellType==1:
	    d2.targetVol = numpy.random.lognormal(1.47,0.62)
	    d2.rad = numpy.random.lognormal(-1.2,0.72)	
	    d2.growthRate =  """ + str(growth_rate) + """	    

    """

    # write script
    f = open("scripts/" + script_name + ".py", "w")
    f.write(script)
    f.close()
