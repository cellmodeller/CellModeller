import sys
import math
import numpy
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
from pyopencl.array import max as device_max
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel
import random
import time


ct_map = {}

class CLBacterium:
    """A rigid body model of bacterial growth implemented using
    OpenCL.
    """

    def __init__(self, simulator,
                 max_substeps=8,
                 max_cells=10000,
                 max_contacts=24,
                 max_planes=1,
                 max_spheres=1,
                 max_sqs=192**2,
                 grid_spacing=5.0,
                 muA=1.0,
                 gamma=10.0,
                 dt=None,
                 cgs_tol=5e-3,
                 jitter_z=True,
                 alternate_divisions=False,
                 printing=True,
                 compNeighbours=False):

        # Should we compute neighbours? (bit slow)
        self.computeNeighbours = compNeighbours

        self.frame_no = 0
        self.simulator = simulator
        self.regulator = None
        
        self.time_begin = time.time()
        self.seconds_elapsed = 0
        self.minutes_elapsed = 0
        self.hours_elapsed = 0

        self.max_cells = max_cells
        self.max_contacts = max_contacts
        self.max_planes = max_planes
        self.max_spheres = max_spheres
        self.max_sqs = max_sqs
        self.grid_spacing = grid_spacing
        self.muA = muA
        self.gamma = gamma
        self.dt = dt
        self.cgs_tol = cgs_tol

        self.max_substeps = max_substeps

        self.n_cells = 0
        self.n_cts = 0
        self.n_planes = 0
        self.n_spheres = 0

        self.next_id = 0

        self.grid_x_min = 0
        self.grid_x_max = 0
        self.grid_y_min = 0
        self.grid_y_max = 0
        self.n_sqs = 0

        self.init_cl()
        #self.init_kernels()
        self.init_data()

        self.parents = {}

        self.jitter_z = jitter_z
        self.alternate_divisions = alternate_divisions
        self.printing = printing
        self.progress_initialised = False
        self.sub_tick_initialised = False

    # Biophysical Model interface
    def reset(self):
        self.n_cells=0
        self.n_cts=0
        self.n_planes=0
        self.n_spheres=0

    def setRegulator(self, regulator):
        self.regulator = regulator
        self.init_kernels()


    def addCell(self, cellState, pos=(0,0,0), dir=(1,0,0), rad=0.5, **kwargs):
        i = cellState.idx
        self.n_cells += 1
        cid = cellState.id
        self.cell_centers[i] = tuple(pos+(0,))
        self.cell_dirs[i] = tuple(dir+(0,))
        self.cell_lens[i] = cellState.length
        self.cell_rads[i] = rad
        self.initCellState(cellState)
        self.set_cells()
        self.calc_cell_geom() # cell needs a volume

    #---
    # Some functions to modify existing cells (e.g. from GUI)
    # Eventually prob better to have a generic editCell() that deals with this stuff
    #
    def moveCell(self, cellState, delta_pos):
        print("cell idx = %d"%cellState.idx)
        i = cellState.idx
        cid = cellState.id
        print("cell center = ")
        print(self.cell_centers[i])
        print("delta_pos")
        print(delta_pos)
        pos = numpy.array(tuple(self.cell_centers[i]))
        pos[0:3] += numpy.array(tuple(delta_pos))
        self.cell_centers[i] = pos
        self.simulator.cellStates[cid].pos = [self.cell_centers[i][j] for j in range(3)]
        self.set_cells()
        self.updateCellState(cellState)

    def addPlane(self, pt, norm, coeff):
        pidx = self.n_planes
        self.n_planes += 1
        self.plane_pts[pidx] = tuple(pt)+(0,)
        self.plane_norms[pidx] = tuple(norm) + (0,)
        self.plane_coeffs[pidx] = coeff
        self.set_planes()

    def addSphere(self, pt, rad, coeff, norm):
        sidx = self.n_spheres
        self.n_spheres += 1
        self.sphere_pts[sidx] = tuple(pt)+(0,)
        self.sphere_rads[sidx] = rad
        self.sphere_coeffs[sidx] = coeff
        self.sphere_norms[sidx] = norm
        self.set_spheres()

    def hasNeighbours(self):
        return False

    def divide(self, parentState, daughter1State, daughter2State, *args, **kwargs):
        self.divide_cell(parentState.idx, daughter1State.idx, daughter2State.idx)
        # Initialise cellState data
        self.initCellState(daughter1State)
        self.initCellState(daughter2State)

    def init_cl(self):
        if self.simulator:
            (self.context, self.queue) = self.simulator.getOpenCL()

    def init_kernels(self):
        """Set up the OpenCL kernels."""
        from pkg_resources import resource_string
        kernel_src = resource_string(__name__, 'CLBacterium.cl').decode()

        self.program = cl.Program(self.context, kernel_src).build(cache_dir=False)
        # Some kernels that seem like they should be built into pyopencl...
        self.vclearf = ElementwiseKernel(self.context, "float8 *v", "v[i]=0.0", "vecclearf")
        self.vcleari = ElementwiseKernel(self.context, "int *v", "v[i]=0", "veccleari")
        self.vadd = ElementwiseKernel(self.context, "float8 *res, const float8 *in1, const float8 *in2",
                                      "res[i] = in1[i] + in2[i]", "vecadd")
        self.vsub = ElementwiseKernel(self.context, "float8 *res, const float8 *in1, const float8 *in2",
                                          "res[i] = in1[i] - in2[i]", "vecsub")
        self.vaddkx = ElementwiseKernel(self.context,
                                            "float8 *res, const float k, const float8 *in1, const float8 *in2",
                                            "res[i] = in1[i] + k*in2[i]", "vecaddkx")
        self.vsubkx = ElementwiseKernel(self.context,
                                            "float8 *res, const float k, const float8 *in1, const float8 *in2",
                                            "res[i] = in1[i] - k*in2[i]", "vecsubkx")
        self.vmulk = ElementwiseKernel(self.context,
                                            "float8 *res, const float k, const float8 *in1",
                                            "res[i] = k*in1[i]", "vecmulk")
        self.vnorm = ElementwiseKernel(self.context,
                                            "float8 *res, const float8 *in1",
                                            "res[i] = dot(in1[i], in1[i]", "vecnorm")


        # cell geometry kernels
        self.calc_cell_area = ElementwiseKernel(self.context, "float* res, float* r, float* l",
                                           "res[i] = 2.f*3.1415927f*r[i]*(2.f*r[i]+l[i])", "cell_area_kern")
        self.calc_cell_vol = ElementwiseKernel(self.context, "float* res, float* r, float* l",
                                          "res[i] = 3.1415927f*r[i]*r[i]*(2.f*r[i]+l[i])", "cell_vol_kern")

        # A dot product as sum of float4 dot products -
        # i.e. like flattening vectors of float8s into big float vectors
        # then computing dot
        # NB. Some openCLs seem not to implement dot(float8,float8) so split
        # into float4's
        self.vdot = ReductionKernel(self.context, numpy.float32, neutral="0",
                reduce_expr="a+b", map_expr="dot(x[i].s0123,y[i].s0123)+dot(x[i].s4567,y[i].s4567)",
                arguments="__global float8 *x, __global float8 *y")
    

    def init_data(self):
        """Set up the data OpenCL will store on the device."""
        # cell data
        cell_geom = (self.max_cells,)
        self.cell_centers = numpy.zeros(cell_geom, vec.float4)
        self.cell_centers_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_dirs = numpy.zeros(cell_geom, vec.float4)
        self.cell_dirs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_lens = numpy.zeros(cell_geom, numpy.float32)
        self.cell_lens_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.pred_cell_centers = numpy.zeros(cell_geom, vec.float4)
        self.pred_cell_centers_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.pred_cell_dirs = numpy.zeros(cell_geom, vec.float4)
        self.pred_cell_dirs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.pred_cell_lens = numpy.zeros(cell_geom, numpy.float32)
        self.pred_cell_lens_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_rads = numpy.zeros(cell_geom, numpy.float32)
        self.cell_rads_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_sqs = numpy.zeros(cell_geom, numpy.int32)
        self.cell_sqs_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)
        self.cell_n_cts = numpy.zeros(cell_geom, numpy.int32)
        self.cell_n_cts_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)
        self.cell_dcenters = numpy.zeros(cell_geom, vec.float4)
        self.cell_dcenters_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_dangs = numpy.zeros(cell_geom, vec.float4)
        self.cell_dangs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_dlens = numpy.zeros(cell_geom, numpy.float32)
        self.cell_dlens_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_target_dlens_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_growth_rates = numpy.zeros(cell_geom, numpy.float32)

        # cell geometry calculated from l and r
        self.cell_areas_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_vols_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_old_vols_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)

        # gridding
        self.sq_inds = numpy.zeros((self.max_sqs,), numpy.int32)
        self.sq_inds_dev = cl_array.zeros(self.queue, (self.max_sqs,), numpy.int32)
        self.sorted_ids = numpy.zeros(cell_geom, numpy.int32)
        self.sorted_ids_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)

        # constraint planes
        plane_geom = (self.max_planes,)
        self.plane_pts = numpy.zeros(plane_geom, vec.float4)
        self.plane_pts_dev = cl_array.zeros(self.queue, plane_geom, vec.float4)
        self.plane_norms = numpy.zeros(plane_geom, vec.float4)
        self.plane_norms_dev = cl_array.zeros(self.queue, plane_geom, vec.float4)
        self.plane_coeffs = numpy.zeros(plane_geom, numpy.float32)
        self.plane_coeffs_dev = cl_array.zeros(self.queue, plane_geom, numpy.float32)

        # constraint spheres
        sphere_geom = (self.max_spheres,)
        self.sphere_pts = numpy.zeros(sphere_geom, vec.float4)
        self.sphere_pts_dev = cl_array.zeros(self.queue, sphere_geom, vec.float4)
        self.sphere_rads = numpy.zeros(sphere_geom, numpy.float32)
        self.sphere_rads_dev = cl_array.zeros(self.queue, sphere_geom, numpy.float32)
        self.sphere_coeffs = numpy.zeros(sphere_geom, numpy.float32)
        self.sphere_coeffs_dev = cl_array.zeros(self.queue, sphere_geom, numpy.float32)
        self.sphere_norms = numpy.zeros(sphere_geom, numpy.float32)
        self.sphere_norms_dev = cl_array.zeros(self.queue, sphere_geom, numpy.float32)

        # contact data
        ct_geom = (self.max_cells, self.max_contacts)
        self.ct_frs = numpy.zeros(ct_geom, numpy.int32)
        self.ct_frs_dev = cl_array.zeros(self.queue, ct_geom, numpy.int32)
        self.ct_tos = numpy.zeros(ct_geom, numpy.int32)
        self.ct_tos_dev = cl_array.zeros(self.queue, ct_geom, numpy.int32)
        self.ct_dists = numpy.zeros(ct_geom, numpy.float32)
        self.ct_dists_dev = cl_array.zeros(self.queue, ct_geom, numpy.float32)
        self.ct_pts = numpy.zeros(ct_geom, vec.float4)
        self.ct_pts_dev = cl_array.zeros(self.queue, ct_geom, vec.float4)
        self.ct_norms = numpy.zeros(ct_geom, vec.float4)
        self.ct_norms_dev = cl_array.zeros(self.queue, ct_geom, vec.float4)
        self.ct_stiff_dev = cl_array.zeros(self.queue, ct_geom, numpy.float32)
        self.ct_overlap_dev = cl_array.zeros(self.queue, ct_geom, numpy.float32)
        self.neighbours = numpy.zeros(ct_geom, numpy.int32)
        self.cell_cts = numpy.zeros(self.max_cells, numpy.int32)

        # where the contacts pointing to this cell are collected
        self.cell_tos = numpy.zeros(ct_geom, numpy.int32)
        self.cell_tos_dev = cl_array.zeros(self.queue, ct_geom, numpy.int32)
        self.n_cell_tos = numpy.zeros(cell_geom, numpy.int32)
        self.n_cell_tos_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)


        # the constructed 'matrix'
        mat_geom = (self.max_cells*self.max_contacts,)
        self.ct_inds = numpy.zeros(mat_geom, numpy.int32)
        self.ct_inds_dev = cl_array.zeros(self.queue, mat_geom, numpy.int32)
        self.ct_reldists = numpy.zeros(mat_geom, numpy.float32)
        self.ct_reldists_dev = cl_array.zeros(self.queue, mat_geom, numpy.float32)

        self.fr_ents = numpy.zeros(mat_geom, vec.float8)
        self.fr_ents_dev = cl_array.zeros(self.queue, mat_geom, vec.float8)
        self.to_ents = numpy.zeros(mat_geom, vec.float8)
        self.to_ents_dev = cl_array.zeros(self.queue, mat_geom, vec.float8)
        

        # vectors and intermediates
        self.deltap = numpy.zeros(cell_geom, vec.float8)
        self.deltap_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
        self.Mx = numpy.zeros(mat_geom, numpy.float32)
        self.Mx_dev = cl_array.zeros(self.queue, mat_geom, numpy.float32)
        self.BTBx = numpy.zeros(cell_geom, vec.float8)
        self.BTBx_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
        self.Minvx_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)

        # CGS intermediates
        self.p_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
        self.Ap_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
        self.res_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
        self.rhs_dev = cl_array.zeros(self.queue, cell_geom, vec.float8)
    

    def load_from_cellstates(self, cell_states):
        for (cid,cs) in list(cell_states.items()):
            i = cs.idx
            self.cell_centers[i] = tuple(cs.pos)+(0,)
            self.cell_dirs[i] = tuple(cs.dir)+(0,)
            self.cell_rads[i] = cs.radius
            self.cell_lens[i] = cs.length
        
        self.n_cells = len(cell_states)
        self.set_cells()
        self.calc_cell_area(self.cell_areas_dev, self.cell_rads_dev, self.cell_lens_dev)
        self.calc_cell_vol(self.cell_vols_dev, self.cell_rads_dev, self.cell_lens_dev)

    def load_test_data(self):
        import CellModeller.Biophysics.BacterialModels.CLData as data
        self.cell_centers.put(list(range(len(data.pos))), data.pos)
        self.cell_dirs.put(list(range(len(data.dirs))), data.dirs)
        self.cell_lens.put(list(range(len(data.lens))), data.lens)
        self.cell_rads.put(list(range(len(data.rads))), data.rads)
        self.n_cells = data.n_cells
        self.set_cells()

    def load_1_cell(self):
        self.cell_centers.put([0], [(0,0,0,0)])
        self.cell_dirs.put([0], [(1,0,0,0)])
        self.cell_lens.put([0], [2.0])
        self.cell_rads.put([0], [0.5])
        self.n_cells = 1
        self.set_cells()


    def load_2_cells(self):
        root2 = numpy.sqrt(2.0)
        self.cell_centers.put([0,1], [(-root2-0.5, 0, 0, 0), (root2+0.5, 0, 0, 0)])
        self.cell_dirs.put([0,1], [(root2/2.0, root2/2.0, 0, 0), (-root2/2.0, root2/2.0, 0, 0)])
        self.cell_lens.put([0,1], [4.0, 4.0])
        self.cell_rads.put([0,1], [0.5, 0.5])
        self.n_cells = 2
        self.set_cells()


    def load_3_cells(self):
        root2 = numpy.sqrt(2.0)
        self.cell_centers.put([0,1,2], [(-root2-0.5, 0, 0, 0), (root2+0.5, 0, 0, 0), (root2+0.5+3.3, 0, 0, 0)])
        self.cell_dirs.put([0,1,2], [(root2/2.0, root2/2.0, 0, 0), (-root2/2.0, root2/2.0, 0, 0), (1, 0, 0, 0)])
        self.cell_lens.put([0,1,2], [3.0, 3.0, 3.0])
        self.cell_rads.put([0,1,2], [0.5, 0.5, 0.5])
        self.n_cells = 3
        self.set_cells()


    def load_3_cells_1_plane(self):
        root2 = numpy.sqrt(2.0)
        self.cell_centers.put([0,1,2], [(-root2-0.5, 0, 0, 0), (root2+0.5, 0, 0, 0), (root2+0.5+3.3, 0, 0, 0)])
        self.cell_dirs.put([0,1,2], [(root2/2.0, root2/2.0, 0, 0), (-root2/2.0, root2/2.0, 0, 0), (1, 0, 0, 0)])
        self.cell_lens.put([0,1,2], [3.0, 3.0, 3.0])
        self.cell_rads.put([0,1,2], [0.5, 0.5, 0.5])
        self.n_cells = 3
        self.set_cells()

        self.n_planes = 1
        self.plane_pts.put([0], [(0, 0, -0.5, 0)])
        self.plane_norms.put([0], [(0, 0, 1, 0)])
        self.plane_coeffs.put([0], [0.5])
        self.set_planes()

    def load_3_cells_2_planes(self):
        root2 = numpy.sqrt(2.0)
        self.cell_centers.put([0,1,2], [(-root2-0.5, 0, 0, 0), (root2+0.5, 0, 0, 0), (root2+0.5+3.3, 0, 0, 0)])
        self.cell_dirs.put([0,1,2], [(root2/2.0, root2/2.0, 0, 0), (-root2/2.0, root2/2.0, 0, 0), (1, 0, 0, 0)])
        self.cell_lens.put([0,1,2], [3.0, 3.0, 3.0])
        self.cell_rads.put([0,1,2], [0.5, 0.5, 0.5])
        self.n_cells = 3
        self.set_cells()

        self.n_planes = 2
        self.plane_pts.put([0,1], [(0, 0, -0.5, 0), (0, 0, 0.5, 0)])
        self.plane_norms.put([0,1], [(0, 0, 1, 0), (0, 0, -1, 0)])
        self.plane_coeffs.put([0,1], [0.5, 0.1])
        self.set_planes()


    def load_1_cell_1_plane(self):
        self.cell_centers.put([0], [(0,0,0,0)])
        self.cell_dirs.put([0], [(1,0,0,0)])
        self.cell_lens.put([0], [3.0])
        self.cell_rads.put([0], [0.5])
        self.n_cells = 1
        self.set_cells()

        self.plane_pts.put([0], [(4, 0, 0, 0)])
        self.plane_norms.put([0], [(-1, 0, 0, 0)])
        self.plane_coeffs.put([0], [0.5])
        self.n_planes = 1
        self.set_planes()


    def load_1024_cells(self):
        d = 32
        for i in range(-d/2,d/2):
            for j in range(-d/2,d/2):
                n = (i+d/2)*d + (j+d/2)
                x = i*3.5 + random.uniform(-0.05,0.05)
                y = j*2.0 + random.uniform(-0.05,0.05)
                th = random.uniform(-0.15, 0.15)
                dir_x = math.cos(th)
                dir_y = math.sin(th)
                self.cell_centers.put([n], [(x, y, 0, 0)])
                self.cell_dirs.put([n], [(dir_x, dir_y, 0, 0)])
                self.cell_lens.put([n], [2])
                self.cell_rads.put([n], 0.5)
        self.n_cells = d*d
        self.set_cells()

    def get_cells(self):
        """Copy cell centers, dirs, lens, and rads from the device."""
        self.cell_centers[0:self.n_cells] = self.cell_centers_dev[0:self.n_cells].get()
        self.cell_dirs[0:self.n_cells] = self.cell_dirs_dev[0:self.n_cells].get()
        self.cell_lens[0:self.n_cells] = self.cell_lens_dev[0:self.n_cells].get()
        self.cell_rads[0:self.n_cells] = self.cell_rads_dev[0:self.n_cells].get()
        self.cell_dlens[0:self.n_cells] = self.cell_dlens_dev[0:self.n_cells].get()
        self.cell_dcenters[0:self.n_cells] = self.cell_dcenters_dev[0:self.n_cells].get()
        self.cell_dangs[0:self.n_cells] = self.cell_dangs_dev[0:self.n_cells].get()

    def set_cells(self):
        """Copy cell centers, dirs, lens, and rads to the device from local."""
        self.cell_centers_dev[0:self.n_cells].set(self.cell_centers[0:self.n_cells])
        self.cell_dirs_dev[0:self.n_cells].set(self.cell_dirs[0:self.n_cells])
        self.cell_lens_dev[0:self.n_cells].set(self.cell_lens[0:self.n_cells])
        self.cell_rads_dev[0:self.n_cells].set(self.cell_rads[0:self.n_cells])
        self.cell_dlens_dev[0:self.n_cells].set(self.cell_dlens[0:self.n_cells])
        self.cell_dcenters_dev[0:self.n_cells].set(self.cell_dcenters[0:self.n_cells])
        self.cell_dangs_dev[0:self.n_cells].set(self.cell_dangs[0:self.n_cells])

    def set_planes(self):
        """Copy plane pts, norms, and coeffs to the device from local."""
        self.plane_pts_dev[0:self.n_planes].set(self.plane_pts[0:self.n_planes])
        self.plane_norms_dev[0:self.n_planes].set(self.plane_norms[0:self.n_planes])
        self.plane_coeffs_dev[0:self.n_planes].set(self.plane_coeffs[0:self.n_planes])


    def set_spheres(self):
        """Copy sphere pts and coeffs to the device from local."""
        self.sphere_pts_dev[0:self.n_spheres].set(self.sphere_pts[0:self.n_spheres])
        self.sphere_rads_dev[0:self.n_spheres].set(self.sphere_rads[0:self.n_spheres])
        self.sphere_coeffs_dev[0:self.n_spheres].set(self.sphere_coeffs[0:self.n_spheres])
        self.sphere_norms_dev[0:self.n_spheres].set(self.sphere_norms[0:self.n_spheres])


    def get_cts(self):
        """Copy contact froms, tos, dists, pts, and norms from the device."""
        self.ct_frs[0:self.n_cts] = self.ct_frs_dev[0:self.n_cts].get()
        self.ct_tos[0:self.n_cts] = self.ct_tos_dev[0:self.n_cts].get()
        self.ct_dists[0:self.n_cts] = self.ct_dists_dev[0:self.n_cts].get()
        self.ct_pts[0:self.n_cts] = self.ct_pts_dev[0:self.n_cts].get()
        self.ct_norms[0:self.n_cts] = self.ct_norms_dev[0:self.n_cts].get()
        self.cell_n_cts[0:self.n_cells] = self.cell_n_cts_dev[0:self.n_cells].get()

    def matrixTest(self):
        x_dev = cl_array.zeros(self.queue, (self.n_cells,), vec.float8)
        Ax_dev = cl_array.zeros(self.queue, (self.n_cells,), vec.float8)
        opstring = ''
        for i in range(self.n_cells):
            x = numpy.zeros((self.n_cells,), vec.float8)
            for j in range(7):
                if j>0:
                    x[i][j-1]=0.0
                x[i][j]=1.0
                x_dev.set(x)
                self.calculate_Ax(Ax_dev, x_dev)
                Ax = Ax_dev.get()
                for ii in range(self.n_cells):
                    for jj in range(7):
                        opstring += str(Ax[ii][jj])
                        if ii!=self.n_cells-1 or jj!=6:
                            opstring = opstring + '\t'
                opstring = opstring + '\n'
        print("MTM")
        print(opstring)
        open('CellModeller/Biophysics/BacterialModels/matrix.mat', 'w').write(opstring)


    def dump_cell_data(self, n):
        import pickle
        filename = 'data/data-%04i.pickle'%n
        outfile = open(filename, 'wb')
        data = (self.n_cells,
                self.cell_centers_dev.get(),
                self.cell_dirs_dev.get(),
                self.cell_lens_dev.get(),
                self.cell_rads_dev.get(),
                self.parents),
        pickle.dump(data, outfile, protocol=-1)

    def dydt(self):
        self.set_cells()

    def finish(self):
        # pull cells from the device and update simulator
        if self.simulator:
            self.get_cells()
            for state in list(self.simulator.cellStates.values()):
                self.updateCellState(state)

    def progress_init(self, dt):
        self.set_cells()
        # NOTE: by default self.dt=None, and time step == simulator time step (dt) 
        if self.dt:
            self.n_ticks = int(math.ceil(dt/self.dt)) 
        else:
            self.n_ticks = 1 
        # print("n_ticks = %d"%(self.n_ticks))
        self.actual_dt = dt / float(self.n_ticks)
        self.progress_initialised = True

    def progress(self):
        if self.n_ticks:
            if self.tick(self.actual_dt):
                self.n_ticks -= 1
            return False
        else:
            return True

    def progress_finalise(self):
        self.frame_no += 1
        self.progress_initialised = False
        self.seconds_elapsed = numpy.float32(time.time() - self.time_begin)
        self.minutes_elapsed = (numpy.float32(self.seconds_elapsed) / 60.0)  
        self.hours_elapsed = (numpy.float32(self.minutes_elapsed) / 60.0)  
        if self.frame_no % 10 == 0:
            print('% 8i    % 8i cells    % 8i contacts    %f hour(s) or %f minute(s) or %f second(s)' % (self.frame_no, self.n_cells, self.n_cts, self.hours_elapsed, self.minutes_elapsed, self.seconds_elapsed))
        # pull cells from the device and update simulator
        if self.simulator:
            self.get_cells()
            # TJR: added incremental construction of this dict to same places as idToIdx - not fully tested
            #idxToId = {idx: id for id, idx in self.simulator.idToIdx.iteritems()}
            # TJR: add flag for this cos a bit time consuming
            if self.computeNeighbours:
                self.updateCellNeighbours(self.simulator.idxToId)
            for state in list(self.simulator.cellStates.values()):
                self.updateCellState(state)

    def step(self, dt):
        """Step forward dt units of time.

        Assumes that:
        cell_centers is up to date when it starts.
        """
        if not self.progress_initialised:
            self.progress_init(dt)
        if self.progress():
            self.progress_finalise()
            return True
        else:
            return False

    def sub_tick_init(self, dt):
        # set target dlens (taken from growth rates set by updateCellStates)
        #self.cell_target_dlens_dev.set(dt*self.cell_growth_rates)
        #self.cell_dlens_dev.set(dt*self.cell_dlens)
        self.cell_dlens_dev[0:self.n_cells].set(dt*self.cell_growth_rates[0:self.n_cells])

        # redefine gridding based on the range of cell positions
        self.cell_centers[0:self.n_cells] = self.cell_centers_dev[0:self.n_cells].get()
        self.update_grid() # we assume local cell_centers is current

        # get each cell into the correct sq and retrieve from the device
        self.bin_cells()

        # sort cells and find sq index starts in the list
        self.cell_sqs = self.cell_sqs_dev[0:self.n_cells].get() # get updated cell sqs
        self.sort_cells()
        self.sorted_ids_dev.set(self.sorted_ids) # push changes to the device
        self.sq_inds_dev.set(self.sq_inds)

        self.n_cts = 0
        self.vcleari(self.cell_n_cts_dev) # clear the accumulated contact count
        self.sub_tick_i=0
        self.sub_tick_initialised=True

    def tick(self, dt):
        if not self.sub_tick_initialised:
            self.sub_tick_init(dt)
        if self.sub_tick(dt):
            self.sub_tick_finalise()
            return True
        else:
            return False

    def sub_tick(self, dt):
        old_n_cts = self.n_cts
        self.predict()
        # find all contacts
        self.find_contacts()
        # place 'backward' contacts in cells
        self.collect_tos()

        self.sub_tick_i += 1
        alpha = 10**(self.sub_tick_i)
        new_cts = self.n_cts - old_n_cts
        if (new_cts>0 or self.sub_tick_i==0) and self.sub_tick_i<self.max_substeps:
            self.build_matrix() # Calculate entries of the matrix
            #print "max cell contacts = %i"%cl_array.max(self.cell_n_cts_dev).get()
            self.CGSSolve(dt, alpha) # invert MTMx to find deltap
            self.add_impulse()
            return False
        else:
            return True

    def sub_tick_finalise(self):
        #print "Substeps = %d"%self.sub_tick_i
        self.integrate()
        self.calc_cell_geom()
        self.sub_tick_initialised=False

    def initCellState(self, state):
        cid = state.id
        i = state.idx
        state.pos = [self.cell_centers[i][j] for j in range(3)]
        state.dir = [self.cell_dirs[i][j] for j in range(3)]
        state.radius = self.cell_rads[i]
        state.length = self.cell_lens[i]
        #for effective growth calulations
        state.oldLen = self.cell_lens[i]
        
        state.volume = state.length # TO DO: do something better here
        pa = numpy.array(state.pos)
        da = numpy.array(state.dir)
        state.ends = (pa-da*state.length*0.5, pa+da*state.length*0.5)
        state.strainRate = 0.0
        #self.cell_dlens[i] = state.growthRate
        state.startVol = state.volume

    def updateCellNeighbours(self, idx2Id):
        ct_tos = self.ct_tos_dev[0:self.n_cells,:].get()
        cell_to_cts = self.cell_n_cts_dev[0:self.n_cells].get()
        cell_cts = numpy.zeros(self.n_cells, numpy.int32)
        for i in range(self.n_cells):
            for j in range(cell_to_cts[i]):
                if ct_tos[i,j]>0:#not a plane contact
                    self.neighbours[i, cell_cts[i]] = idx2Id[ct_tos[i,j]]
                    cell_cts[i] += 1
                    self.neighbours[ct_tos[i,j], cell_cts[ct_tos[i,j]]] = idx2Id[i]
                    cell_cts[ct_tos[i,j]] += 1
        self.cell_cts = cell_cts

    def updateCellState(self, state):
        cid = state.id
        i = state.idx

        state.vel = [self.cell_centers[i][j]-state.pos[j] for j in range(3)]
        state.pos = [self.cell_centers[i][j] for j in range(3)]
        state.dir = [self.cell_dirs[i][j] for j in range(3)]
        state.radius = self.cell_rads[i]
        state.length = self.cell_lens[i]
        state.strainRate = (state.length - state.oldLen)/state.oldLen

        #currently the effective growth rate is calculated over the entire history of the cell
        state.effGrowth = state.effGrowth * state.cellAge + state.strainRate * state.oldLen
        state.cellAge += 1
        state.effGrowth = state.effGrowth / state.cellAge
        state.oldLen = state.length

        state.neighbours = [] #clear contacts
        if self.computeNeighbours: #populate cellstate.neighbours
            for n in range(self.cell_cts[i]):
                if self.neighbours[i,n] not in state.neighbours:
                    state.neighbours.append(self.neighbours[i,n]) #ids of all cells in physical contact
        state.cts = len(state.neighbours)

        state.volume = state.length # TO DO: do something better here
        pa = numpy.array(state.pos)
        da = numpy.array(state.dir)
        state.ends = (pa-da*state.length*0.5, pa+da*state.length*0.5)
        # Length vel is linearisation of exponential growth
        self.cell_growth_rates[i] = state.growthRate*state.length


    def update_grid(self):
        """Update our grid_(x,y)_min, grid_(x,y)_max, and n_sqs.

        Assumes that our copy of cell_centers is current.
        """
        coords = self.cell_centers.view(numpy.float32).reshape((self.max_cells, 4))

        x_coords = coords[:,0]
        min_x_coord = x_coords.min()
        max_x_coord = x_coords.max()
        self.grid_x_min = int(math.floor(min_x_coord / self.grid_spacing))
        self.grid_x_max = int(math.ceil(max_x_coord / self.grid_spacing))
        if self.grid_x_min == self.grid_x_max:
            self.grid_x_max += 1

        y_coords = coords[:,1]
        min_y_coord = y_coords.min()
        max_y_coord = y_coords.max()
        self.grid_y_min = int(math.floor(min_y_coord / self.grid_spacing))
        self.grid_y_max = int(math.ceil(max_y_coord / self.grid_spacing))
        if self.grid_y_min == self.grid_y_max:
            self.grid_y_max += 1

        self.n_sqs = (self.grid_x_max-self.grid_x_min)*(self.grid_y_max-self.grid_y_min)


    def bin_cells(self):
        """Call the bin_cells kernel.

        Assumes cell_centers is current on the device.

        Calculates cell_sqs.
        """
        self.program.bin_cells(self.queue,
                               (self.n_cells,),
                               None,
                               numpy.int32(self.grid_x_min),
                               numpy.int32(self.grid_x_max),
                               numpy.int32(self.grid_y_min),
                               numpy.int32(self.grid_y_max),
                               numpy.float32(self.grid_spacing),
                               self.cell_centers_dev.data,
                               self.cell_sqs_dev.data).wait()


    def sort_cells(self):
        """Sort the cells by grid square and find the start of each
        grid square's cells in that list.

        Assumes that the local copy of cell_sqs is current.

        Calculates local sorted_ids and sq_inds.
        """
        self.sorted_ids.put(numpy.arange(self.n_cells), numpy.argsort(self.cell_sqs[:self.n_cells]))
        self.sorted_ids_dev[0:self.n_cells].set(self.sorted_ids[0:self.n_cells])

        # find the start of each sq in the list of sorted cell ids and send to the device
        sorted_sqs = numpy.sort(self.cell_sqs[:self.n_cells])
        self.sq_inds.put(numpy.arange(self.n_sqs), numpy.searchsorted(sorted_sqs, numpy.arange(self.n_sqs), side='left'))
        self.sq_inds_dev.set(self.sq_inds)


    def find_contacts(self, predict=True):
        """Call the find_contacts kernel.

        Assumes that cell_centers, cell_dirs, cell_lens, cell_rads,
        cell_sqs, cell_dcenters, cell_dlens, cell_dangs,
        sorted_ids, and sq_inds are current on the device.

        Calculates cell_n_cts, ct_frs, ct_tos, ct_dists, ct_pts,
        ct_norms, ct_reldists, and n_cts.
        """
        if predict:
            centers = self.pred_cell_centers_dev
            dirs = self.pred_cell_dirs_dev
            lens = self.pred_cell_lens_dev
        else:
            centers = self.cell_centers_dev
            dirs = self.cell_dirs_dev
            lens = self.cell_lens_dev

        self.program.find_plane_contacts(self.queue,
                                         (self.n_cells,),
                                         None,
                                         numpy.int32(self.max_cells),
                                         numpy.int32(self.max_contacts),
                                         numpy.int32(self.n_planes),
                                         self.plane_pts_dev.data,
                                         self.plane_norms_dev.data,
                                         self.plane_coeffs_dev.data,
                                         centers.data,
                                         dirs.data,
                                         lens.data,
                                         self.cell_rads_dev.data,
                                         self.cell_n_cts_dev.data,
                                         self.ct_frs_dev.data,
                                         self.ct_tos_dev.data,
                                         self.ct_dists_dev.data,
                                         self.ct_pts_dev.data,
                                         self.ct_norms_dev.data,
                                         self.ct_reldists_dev.data,
                                         self.ct_stiff_dev.data).wait()

        self.program.find_sphere_contacts(self.queue,
                                         (self.n_cells,),
                                         None,
                                         numpy.int32(self.max_cells),
                                         numpy.int32(self.max_contacts),
                                         numpy.int32(self.n_spheres),
                                         self.sphere_pts_dev.data,
                                         self.sphere_coeffs_dev.data,
                                         self.sphere_rads_dev.data,
                                         self.sphere_norms_dev.data,
                                         centers.data,
                                         dirs.data,
                                         lens.data,
                                         self.cell_rads_dev.data,
                                         self.cell_n_cts_dev.data,
                                         self.ct_frs_dev.data,
                                         self.ct_tos_dev.data,
                                         self.ct_dists_dev.data,
                                         self.ct_pts_dev.data,
                                         self.ct_norms_dev.data,
                                         self.ct_reldists_dev.data,
                                         self.ct_stiff_dev.data).wait()

        self.program.find_contacts(self.queue,
                                   (self.n_cells,),
                                   None,
                                   numpy.int32(self.max_cells),
                                   numpy.int32(self.n_cells),
                                   numpy.int32(self.grid_x_min),
                                   numpy.int32(self.grid_x_max),
                                   numpy.int32(self.grid_y_min),
                                   numpy.int32(self.grid_y_max),
                                   numpy.int32(self.n_sqs),
                                   numpy.int32(self.max_contacts),
                                   centers.data,
                                   dirs.data,
                                   lens.data,
                                   self.cell_rads_dev.data,
                                   self.cell_sqs_dev.data,
                                   self.sorted_ids_dev.data,
                                   self.sq_inds_dev.data,
                                   self.cell_n_cts_dev.data,
                                   self.ct_frs_dev.data,
                                   self.ct_tos_dev.data,
                                   self.ct_dists_dev.data,
                                   self.ct_pts_dev.data,
                                   self.ct_norms_dev.data,
                                   self.ct_reldists_dev.data,
                                   self.ct_stiff_dev.data,
                                   self.ct_overlap_dev.data).wait()

        # set dtype to int32 so we don't overflow the int32 when summing
        #self.n_cts = self.cell_n_cts_dev.get().sum(dtype=numpy.int32)
        self.n_cts = cl_array.sum(self.cell_n_cts_dev[0:self.n_cells]).get()


    def collect_tos(self):
        """Call the collect_tos kernel.

        Assumes that cell_sqs, sorted_ids, sq_inds, cell_n_cts,
        ct_frs, and ct_tos are current on the device.

        Calculates cell_tos and n_cell_tos.
        """
        self.program.collect_tos(self.queue,
                                 (self.n_cells,),
                                 None,
                                 numpy.int32(self.max_cells),
                                 numpy.int32(self.n_cells),
                                 numpy.int32(self.grid_x_min),
                                 numpy.int32(self.grid_x_max),
                                 numpy.int32(self.grid_y_min),
                                 numpy.int32(self.grid_y_max),
                                 numpy.int32(self.n_sqs),
                                 numpy.int32(self.max_contacts),
                                 self.cell_sqs_dev.data,
                                 self.sorted_ids_dev.data,
                                 self.sq_inds_dev.data,
                                 self.cell_n_cts_dev.data,
                                 self.ct_frs_dev.data,
                                 self.ct_tos_dev.data,
                                 self.cell_tos_dev.data,
                                 self.n_cell_tos_dev.data).wait()


    def build_matrix(self):
        """Build the matrix so we can calculate M^TMx = Ax.

        Assumes cell_centers, cell_dirs, cell_lens, cell_rads,
        ct_inds, ct_frs, ct_tos, ct_dists, and ct_norms are current on
        the device.

        Calculates fr_ents and to_ents.
        """
        self.program.build_matrix(self.queue,
                                  (self.n_cells, self.max_contacts),
                                  None,
                                  numpy.int32(self.max_contacts),
                                  self.pred_cell_centers_dev.data,
                                  self.pred_cell_dirs_dev.data,
                                  self.pred_cell_lens_dev.data,
                                  self.cell_rads_dev.data,
                                  self.cell_n_cts_dev.data,
                                  self.ct_frs_dev.data,
                                  self.ct_tos_dev.data,
                                  self.ct_pts_dev.data,
                                  self.ct_norms_dev.data,
                                  self.fr_ents_dev.data,
                                  self.to_ents_dev.data,
                                  self.ct_stiff_dev.data).wait()
    

    def calculate_Ax(self, Ax, x, dt, alpha):

        self.program.calculate_Bx(self.queue,
                                  (self.n_cells, self.max_contacts),
                                  None,
                                  numpy.int32(self.max_contacts),
                                  self.ct_frs_dev.data,
                                  self.ct_tos_dev.data,
                                  self.fr_ents_dev.data,
                                  self.to_ents_dev.data,
                                  x.data,
                                  self.Mx_dev.data).wait()
        self.program.calculate_BTBx(self.queue,
                                    (self.n_cells,),
                                    None,
                                    numpy.int32(self.max_contacts),
                                    self.cell_n_cts_dev.data,
                                    self.n_cell_tos_dev.data,
                                    self.cell_tos_dev.data,
                                    self.fr_ents_dev.data,
                                    self.to_ents_dev.data,
                                    self.Mx_dev.data,
                                    Ax.data).wait()
        # Tikhonov test
        #self.vaddkx(Ax, numpy.float32(0.01), Ax, x)

        # Energy mimizing regularization
        self.program.calculate_Mx(self.queue,
                                      (self.n_cells,),
                                      None,
                                      numpy.float32(self.muA),
                                      numpy.float32(self.gamma),
                                      self.cell_dirs_dev.data,
                                      self.cell_lens_dev.data,
                                      self.cell_rads_dev.data,
                                      x.data,
                                      self.Mx_dev.data).wait()

        #this was altered from dt*reg_param
        #self.vaddkx(Ax, self.gamma, Ax, self.Mx_dev).wait()
        #self.vaddkx(Ax, alpha, self.Mx_dev, Ax).wait()
        self.vaddkx(Ax, 1/self.gamma, Ax, self.Mx_dev).wait()
        # 1/math.sqrt(self.n_cells) removed from the reg_param NB
    
        #print(self.Minvx_dev)

    def CGSSolve(self, dt, alpha, substep=False):
        # Solve A^TA\deltap=A^Tb (Ax=b)

        # There must be a way to do this using built in pyopencl - what
        # is it?!
        self.vclearf(self.deltap_dev[0:self.n_cells])
        self.vclearf(self.rhs_dev[0:self.n_cells])

        # put M^T n^Tv_rel in rhs (b)
        self.program.calculate_BTBx(self.queue,
                                    (self.n_cells,),
                                    None,
                                    numpy.int32(self.max_contacts),
                                    self.cell_n_cts_dev.data,
                                    self.n_cell_tos_dev.data,
                                    self.cell_tos_dev.data,
                                    self.fr_ents_dev.data,
                                    self.to_ents_dev.data,
                                    self.ct_reldists_dev.data,
                                    self.rhs_dev.data).wait()


        # res = b-Ax
        self.calculate_Ax(self.BTBx_dev, self.deltap_dev, dt, alpha)
        self.vsub(self.res_dev[0:self.n_cells], self.rhs_dev[0:self.n_cells], self.BTBx_dev[0:self.n_cells])

        # p = res
        cl.enqueue_copy(self.queue, self.p_dev[0:self.n_cells].data, self.res_dev[0:self.n_cells].data)

        # rsold = l2norm(res)
        rsold = self.vdot(self.res_dev[0:self.n_cells], self.res_dev[0:self.n_cells]).get()
        rsfirst = rsold
        if math.sqrt(rsold/self.n_cells) < self.cgs_tol:
            if self.printing and self.frame_no%10==0:
                print('% 5i'%self.frame_no + '% 6i cells  % 6i cts  % 6i iterations  residual = %f' % (self.n_cells,
            self.n_cts, 0, math.sqrt(rsold/self.n_cells)))
            return (0.0, rsold)

        # iterate
        # max iters = matrix dimension = 7 (dofs) * num cells
        #dying=False
        max_iters = self.n_cells*7
            
        for iter in range(max_iters):
            # Ap
            self.calculate_Ax(self.Ap_dev[0:self.n_cells], self.p_dev[0:self.n_cells], dt, alpha)

            # p^TAp
            pAp = self.vdot(self.p_dev[0:self.n_cells], self.Ap_dev[0:self.n_cells]).get()

            # alpha = rsold/p^TAp
            alpha = numpy.float32(rsold/pAp)

            # x = x + alpha*p, x=self.disp
            self.vaddkx(self.deltap_dev[0:self.n_cells], alpha, self.deltap_dev[0:self.n_cells], self.p_dev[0:self.n_cells])

            # res = res - alpha*Ap
            self.vsubkx(self.res_dev[0:self.n_cells], alpha, self.res_dev[0:self.n_cells], self.Ap_dev[0:self.n_cells])

            # rsnew = l2norm(res)
            rsnew = self.vdot(self.res_dev[0:self.n_cells], self.res_dev[0:self.n_cells]).get()

            # Test for convergence
            if math.sqrt(rsnew/self.n_cells) < self.cgs_tol:
            #if math.sqrt(rsnew/rsfirst) < self.cgs_tol:
                break

            # Stopped converging -> terminate
            #if rsnew/rsold>2.0:
            #    break

            # p = res + rsnew/rsold *p
            self.vaddkx(self.p_dev[0:self.n_cells], numpy.float32(rsnew/rsold), self.res_dev[0:self.n_cells], self.p_dev[0:self.n_cells])

            rsold = rsnew
            #print '        ',iter,rsold

        if self.printing and self.frame_no%10==0:
            print('% 5i'%self.frame_no + '% 6i cells  % 6i cts  % 6i iterations  residual = %f' % (self.n_cells, self.n_cts, iter+1, math.sqrt(rsnew/self.n_cells)))
        return (iter+1, math.sqrt(rsnew/self.n_cells))


    def predict(self):
        """Predict cell centers, dirs, lens for a timestep dt based
        on the current velocities.

        Assumes cell_centers, cell_dirs, cell_lens, cell_rads, and
        cell_dcenters, cell_dangs, cell_dlens are current on the device.

        Calculates new pred_cell_centers, pred_cell_dirs, pred_cell_lens.
        """
        self.program.predict(self.queue,
                             (self.n_cells,),
                             None,
                             self.cell_centers_dev.data,
                             self.cell_dirs_dev.data,
                             self.cell_lens_dev.data,
                             self.cell_dcenters_dev.data,
                             self.cell_dangs_dev.data,
                             self.cell_dlens_dev.data,
                             self.pred_cell_centers_dev.data,
                             self.pred_cell_dirs_dev.data,
                             self.pred_cell_lens_dev.data).wait()

    def integrate(self):
        """Integrates cell centers, dirs, lens for a timestep dt based
        on the current deltap.

        Assumes cell_centers, cell_dirs, cell_lens, cell_rads, and
        deltap are current on the device.

        Calculates new cell_centers, cell_dirs, cell_lens.
        """
        self.program.integrate(self.queue,
                               (self.n_cells,),
                               None,
                               self.cell_centers_dev.data,
                               self.cell_dirs_dev.data,
                               self.cell_lens_dev.data,
                               self.cell_dcenters_dev.data,
                               self.cell_dangs_dev.data,
                               self.cell_dlens_dev.data).wait()

    def add_impulse(self):
        self.program.add_impulse(self.queue, (self.n_cells,), None,
                                 numpy.float32(self.muA),
                                 numpy.float32(self.gamma),
                                 self.deltap_dev.data,
                                 self.cell_dirs_dev.data,
                                 self.cell_lens_dev.data,
                                 self.cell_rads_dev.data,
                                 self.cell_dcenters_dev.data,
                                 self.cell_dangs_dev.data,
                                 self.cell_target_dlens_dev.data,
                                 self.cell_dlens_dev.data).wait()

    def divide_cell(self, i, d1i, d2i):
        """Divide a cell into two equal sized daughter cells.

        Fails silently if we're out of cells.

        Assumes our local copy of cells is current.

        Calculates new cell_centers, cell_dirs, cell_lens, and cell_rads.
        """
        if self.n_cells >= self.max_cells:
            return
        # idxs of the two new cells
        a = d1i
        b = d2i

        # seems to be making shallow copies without the tuple calls
        parent_center = tuple(self.cell_centers[i])
        parent_dir = tuple(self.cell_dirs[i])
        parent_rad = self.cell_rads[i]
        parent_len = self.cell_lens[i]

        daughter_len = parent_len/2.0 - parent_rad #- 0.025
        daughter_offset = daughter_len/2.0 + parent_rad
        center_offset = tuple([parent_dir[k]*daughter_offset for k in range(4)])

        self.cell_centers[a] = tuple([(parent_center[k] - center_offset[k]) for k in range(4)])
        self.cell_centers[b] = tuple([(parent_center[k] + center_offset[k]) for k in range(4)])

        if not self.alternate_divisions:
            cdir = numpy.array(parent_dir)
            jitter = numpy.random.uniform(-0.001,0.001,3)
            if not self.jitter_z: jitter[2] = 0.0
            cdir[0:3] += jitter
            cdir /= numpy.linalg.norm(cdir)
            self.cell_dirs[a][0] = cdir[0]
            self.cell_dirs[a][1] = cdir[1]
            self.cell_dirs[a][2] = cdir[2]

            cdir = numpy.array(parent_dir)
            jitter = numpy.random.uniform(-0.001,0.001,3)
            if not self.jitter_z: jitter[2] = 0.0
            cdir[0:3] += jitter
            cdir /= numpy.linalg.norm(cdir)
            self.cell_dirs[b][0] = cdir[0]
            self.cell_dirs[b][1] = cdir[1]
            self.cell_dirs[b][2] = cdir[2]
        else:
            cdir = numpy.array(parent_dir)
            tmp = cdir[0]
            cdir[0] = -cdir[1]
            cdir[1] = tmp
            self.cell_dirs[a] = cdir
            self.cell_dirs[b] = cdir


        self.cell_lens[a] = daughter_len
        self.cell_lens[b] = daughter_len
        self.cell_rads[a] = parent_rad
        self.cell_rads[b] = parent_rad

        self.n_cells += 1

        self.parents[b] = a

        vols = self.cell_vols_dev[0:self.n_cells].get()
        daughter_vol = vols[i] / 2.0
        vols[a] = daughter_vol
        vols[b] = daughter_vol
        self.cell_vols_dev[0:self.n_cells].set(vols)

        # Inherit velocities from parent (conserve momentum)
        parent_dlin = self.cell_dcenters[i]
        self.cell_dcenters[a] = parent_dlin
        self.cell_dcenters[b] = parent_dlin
        parent_dang = self.cell_dangs[i]
        self.cell_dangs[a] = parent_dang
        self.cell_dangs[b] = parent_dang


        #return indices of daughter cells
        return (a,b)


    def calc_cell_geom(self):
        """Calculate cell geometry using lens/rads on card."""
        # swap cell vols and cell_vols old
        tmp = self.cell_old_vols_dev[0:self.n_cells]
        self.cell_old_vols_dev[0:self.n_cells] = self.cell_vols_dev[0:self.n_cells]
        self.cell_vols_dev[0:self.n_cells] = tmp
        # update geometry
        self.calc_cell_area(self.cell_areas_dev[0:self.n_cells], \
                            self.cell_rads_dev[0:self.n_cells], \
                            self.cell_lens_dev[0:self.n_cells])
        self.calc_cell_vol(self.cell_vols_dev[0:self.n_cells], \
                            self.cell_rads_dev[0:self.n_cells], \
                            self.cell_lens_dev[0:self.n_cells])


    def profileGrid(self):
        if self.n_cts==0:
            return
        import time
        t1 = time.clock()
        for i in range(1000):
            # redefine gridding based on the range of cell positions
            self.cell_centers = self.cell_centers_dev.get()
            self.update_grid() # we assume local cell_centers is current

            # get each cell into the correct sq and retrieve from the device
            self.bin_cells()

            # sort cells and find sq index starts in the list
            self.cell_sqs = self.cell_sqs_dev.get() # get updated cell sqs
            self.sort_cells()
            self.sorted_ids_dev.set(self.sorted_ids) # push changes to the device
            self.sq_inds_dev.set(self.sq_inds)
        t2 = time.clock()
        print("Grid stuff timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
        open("grid_prof","a").write( "%i, %i, %f\n"%(self.n_cells,self.n_cts,(t2-t1)*0.001) )


    def profileFindCts(self):
        if self.n_cts==0:
            return
        import time
        t1 = time.clock()
        dt = 0.005
        for i in range(1000):
            self.n_cts = 0
            self.vcleari(self.cell_n_cts_dev) # clear the accumulated contact count
            self.predict(dt)
            # find all contacts
            self.find_contacts(dt)
            # place 'backward' contacts in cells
            self.collect_tos()

            # compact the contacts so we can dispatch only enough threads
            # to deal with each
            self.ct_frs = self.ct_frs_dev.get()
            self.ct_tos = self.ct_tos_dev.get()
            self.compact_cts()
            self.ct_inds_dev.set(self.ct_inds)
        t2 = time.clock()
        print("Find contacts timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
        open("findcts_prof","a").write( "%i, %i, %f\n"%(self.n_cells,self.n_cts,(t2-t1)*0.001) )

    def profileCGS(self):
        if self.n_cts==0:
            return
        import time
        t1 = time.clock()
        dt = 0.005
        for i in range(1000):
            self.build_matrix(dt) # Calculate entries of the matrix
            (iters, res) = self.CGSSolve()
            print("cgs prof: iters=%i, res=%f"%(iters,res))
        t2 = time.clock()
        print("CGS timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
        open("cgs_prof","a").write( "%i, %i, %i, %f\n"%(self.n_cells,self.n_cts,iters,(t2-t1)*0.001) )


