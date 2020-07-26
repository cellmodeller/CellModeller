from fenics import *
import numpy as np
from dolfin import *
from pyopencl.array import vec

class FEMDiffusion:
    """
    Diffusion solver using Fenics finite element method.
    """
    def __init__(self, sim, mesh_file, pvd_file, nSignals, diffusion_rate, dt, regul=None):
        self.isGrid = False
        self.nSignals = nSignals

        if pvd_file:
            self.file = File(pvd_file)
        else:
            self.file = None

        self.mesh = Mesh(mesh_file)
        self.dV = 1.
        self.dt = dt
        self.diffusion_rate = diffusion_rate

        # Build FEM system
        self.build_system()
        # Initial value of solution
        self.u = Function(self.V)

        self.regul = regul 
        self.setCellStates(sim.cellStates)

    def build_system(self):
        # Set up FEM problem and function space
        self.V = FunctionSpace(self.mesh, "CG", 1)
        self.V_vec = VectorFunctionSpace(self.mesh, "CG", 1)
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        u0 = Constant(0)
        f = Constant(0)
        self.t = 0
        a = u * v * dx + self.dt * self.diffusion_rate * inner(grad(u), grad(v)) * dx 
        L = u0*v*dx + self.dt*f*v*dx
        self.bc = DirichletBC(self.V, Constant(0), DomainBoundary())
        self.A, self.b = assemble_system(a, L, self.bc)

    def saveData(self, data):
        sig_data = {
                'u': self.u,
                'mesh': self.mesh
                }
        data.update(sig_data)
        return data

    def loadData(self, data):
        self.mesh = data['mesh']
        self.build_system()
        self.u = data['u']

    def setCellStates(self, cs):
        self.cellStates = cs

    def setBiophysics(self, biophysics):
        self.biophys = biophysics

    def setRegulator(self, regul):
        self.regul = regul
    
    def addCell(self, cellState):
        pass

    def gradient(self, gx_dev, gy_dev, gz_dev):
        g = project(grad(self.u),self.V_vec)
        gradx = [g(c.pos[0], c.pos[1], c.pos[2])[0] for id,c in self.cellStates.items()]
        grady = [g(c.pos[0], c.pos[1], c.pos[2])[1] for id,c in self.cellStates.items()]
        gradz = [g(c.pos[0], c.pos[1], c.pos[2])[2] for id,c in self.cellStates.items()]
        gx_dev.set(np.array(gradx, dtype=np.float32))
        gy_dev.set(np.array(grady, dtype=np.float32))
        gz_dev.set(np.array(gradz, dtype=np.float32))

    def add_point_source(self, pos, rate):
        x, y, z = pos
        delta = PointSource(self.V, Point(x, y, z), rate)
        delta.apply(self.b)

    def add_point_sources(self, sources):
        source = PointSource(self.V, sources)
        source.apply(self.b)

    def signals(self, cellSigLevels_dev):
        sigs = [[self.u(c.pos[0], c.pos[1], c.pos[2])] for id,c in self.cellStates.items()]
        sigs = np.array(sigs, dtype=np.float32)
        cellSigLevels_dev.set(sigs)

    def initSignalLevels(self, levels):
        pass

    def step(self, dt):
        # Solve FEM system
        solve(self.A, self.u.vector(), self.b)
        # Reset point sources
        self.b.zero()
        self.bc.apply(self.b)
        self.t += dt
        if self.file:
            self.file << (self.u, self.t)
