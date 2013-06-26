from CellModeller.Biophysics.BacterialModels import CGSBacteriumExt

class CGSBacterium:
    def __init__(self, sim):
	self.solver = CGSBacteriumExt.Solver()
	self.sim = sim

    def setRegulator(self, reg):
	self.reg = reg

    def setSignalling(self, sig):
	self.signalling = sig

    def updateCellState(self, cellState):
        cid = cellState.id
        cellState.pos = self.solver.pos(cid)
	cellState.volume = self.solver.volume(cid)
        cellState.length = self.solver.length(cid)
        cellState.radius = self.solver.radius(cid)
        cellState.ends = self.solver.ends(cid)
        self.solver.setGrowthRate(cid, cellState.growthRate)

    def initCellState(self, cellState):
        cid = cellState.id
        cellState.pos = self.solver.pos(cid)
	cellState.volume = self.solver.volume(cid)
        cellState.length = self.solver.length(cid)
        cellState.radius = self.solver.radius(cid)
        cellState.ends = self.solver.ends(cid)
        self.solver.setGrowthRate(cid, cellState.growthRate)
	cellState.startVol = cellState.volume

    def addCell(self, cellState):
        # puts a bug at (0,0,0)
        cid = cellState.id
        self.solver.addCell(cid)
        # initialise geometry data
        self.initCellState(cellState)
	
    def divide(self, state, d1State, d2State, *args, **kw):
        f1 = 1
        f2 = 1
        if kw.has_key('f1') and kw.has_key('f2'):
            f1 = kw['f1']
            f2 = kw['f2']
        self.solver.divideCell(state.id, d1State.id, d2State.id, f1, f2)
        self.initCellState(d1State)
        self.initCellState(d2State)
                
    def step(self, dt):
	self.solver.step(dt)
	for (cid,c) in self.sim.cellStates.items():
            self.updateCellState(c)
