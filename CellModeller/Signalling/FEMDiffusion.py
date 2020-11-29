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
        self.u = []
        for sig in range(self.nSignals):
            self.u.append(Function(self.V[sig]))

        self.regul = regul 
        self.setCellStates(sim.cellStates)

    def build_system(self):
        self.V = []
        self.V_vec = []
        self.bc = []
        self.A = []
        self.b = []
        self.t = 0
        for s in range(self.nSignals):
            # Set up FEM problem and function space
            V = FunctionSpace(self.mesh, "CG", 1)
            V_vec = VectorFunctionSpace(self.mesh, "CG", 1)
            self.V.append(V)
            self.V_vec.append(V_vec)
            u = TrialFunction(V)
            v = TestFunction(V)
            u0 = Constant(0)
            f = Constant(0)
            a = u * v * dx + self.dt * self.diffusion_rate[s] * inner(grad(u), grad(v)) * dx 
            L = u0*v*dx + self.dt*f*v*dx
            bc = DirichletBC(V, Constant(0), DomainBoundary())
            self.bc.append(bc)
            A, b = assemble_system(a, L, bc)
            self.A.append(A)
            self.b.append(b)

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
        gx = []
        gy = []
        gz = []
        for signal in range(self.nSignals):
            g = project(grad(self.u[signal]),self.V_vec[signal])
            gradx = [g(c.pos[0], c.pos[1], c.pos[2])[0] for id,c in self.cellStates.items()]
            grady = [g(c.pos[0], c.pos[1], c.pos[2])[1] for id,c in self.cellStates.items()]
            gradz = [g(c.pos[0], c.pos[1], c.pos[2])[2] for id,c in self.cellStates.items()]
            gx.append(gradx)
            gy.append(grady)
            gz.append(gradz)
        gx_dev.set(np.array(gx, dtype=np.float32))
        gy_dev.set(np.array(gy, dtype=np.float32))
        gz_dev.set(np.array(gz, dtype=np.float32))

    def add_point_source(self, pos, rate, signal):
        x, y, z = pos
        delta = PointSource(self.V[signal], Point(x, y, z), rate)
        delta.apply(self.b[signal])

    def add_point_sources(self, sources, signal):
        source = PointSource(self.V[signal], sources)
        source.apply(self.b[signal])

    def signals(self, cellSigLevels_dev):
        sigs = [[u(c.pos[0], c.pos[1], c.pos[2]) for u in self.u] for id,c in self.cellStates.items()]
        sigs = np.array(sigs, dtype=np.float32)
        cellSigLevels_dev.set(sigs)

    def initSignalLevels(self, levels):
        pass

    def step(self, dt):
        self.t += dt
        for sig in range(self.nSignals):
            # Solve FEM system
            solve(self.A[sig], self.u[sig].vector(), self.b[sig])
            # Reset point sources
            self.b[sig].zero()
            self.bc[sig].apply(self.b[sig])
            if self.file:
                self.file << (self.u[sig], self.t)
