import sys
import math
import numpy
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
from pyopencl.array import max as device_max
from pyopencl.elementwise import ElementwiseKernel
from pyopencl.reduction import ReductionKernel
from pyopencl.clmath import cos, sin
import random
import time
from pyopencl.clrandom import PhiloxGenerator


ct_map = {}

class CLSPP:
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
                 gamma_s=1.0,
                 Fm=1.,
                 Ws=1.,
                 Wc=0,
                 fcil=0,
                 ftax=0,
                 forg=0,
                 D=1.,
                 dt=None,
                 cgs_tol=5e-3,
                 jitter_z=True,
                 alternate_divisions=False,
                 printing=True,
                 compNeighbours=False,
                 spherical=True,
                 sphere_radius=50):

        # Is the simulaiton on a sphere?
        self.spherical = spherical
        self.sphere_radius = sphere_radius

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
        self.gamma_s = gamma_s
        self.Fm = Fm
        self.Ws = Ws
        self.Wc = Wc
        self.fcil = fcil
        self.ftax = ftax
        self.forg = forg
        self.D = D
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

        # Random number generator
        self.rand = PhiloxGenerator(self.context)

    def __del__(self):
        self.cell_centers_dev.data.release()
        self.cell_angs_dev.data.release()
        self.cell_dirs_dev.data.release()
        self.pred_cell_centers_dev.data.release()
        self.pred_cell_dirs_dev.data.release()
        self.cell_rads_dev.data.release()
        self.cell_sqs_dev.data.release()
        self.cell_n_cts_dev.data.release()
        self.cell_dcenters_dev.data.release()
        self.avg_neighbour_dir_dev.data.release()

        self.cell_areas_dev.data.release()
        self.cell_vols_dev.data.release()
        # gridding
        self.sq_inds_dev.data.release()
        self.sorted_ids_dev.data.release()

        # constraint planes
        self.plane_pts_dev.data.release()
        self.plane_norms_dev.data.release()
        self.plane_coeffs_dev.data.release()

        # constraint spheres
        self.sphere_pts_dev.data.release()
        self.sphere_rads_dev.data.release()
        self.sphere_coeffs_dev.data.release()
        self.sphere_norms_dev.data.release()

        # contact data
        self.ct_frs_dev.data.release()
        self.ct_tos_dev.data.release()
        self.ct_dists_dev.data.release()
        self.ct_pts_dev.data.release()
        self.ct_norms_dev.data.release()
        self.ct_stiff_dev.data.release()
        self.ct_overlap_dev.data.release()

        # where the contacts pointing to this cell are collected
        self.cell_tos_dev.data.release()
        self.n_cell_tos_dev.data.release()


        # the constructed 'matrix'
        self.ct_inds_dev.data.release()
        self.ct_reldists_dev.data.release()

        self.fr_ents_dev.data.release()
        self.to_ents_dev.data.release()
        

        # vectors and intermediates
        self.deltap_dev.data.release()
        self.Mx_dev.data.release()
        self.BTBx_dev.data.release()
        self.Minvx_dev.data.release()

        # CGS intermediates
        self.p_dev.data.release()
        self.Ap_dev.data.release()
        self.res_dev.data.release()
        self.rhs_dev.data.release()

    # Biophysical Model interface
    def reset(self):
        self.n_cells=0
        self.n_cts=0
        self.n_planes=0
        self.n_spheres=0

    def setRegulator(self, regulator):
        self.regulator = regulator
        self.init_kernels()


    def addCell(self, cellState, pos=(0,0,0), dir=(1,0,0), rad=1., **kwargs):
        i = cellState.idx
        self.n_cells += 1
        cid = cellState.id
        self.cell_centers[i] = tuple(pos+(0,))
        self.cell_dirs[i] = tuple(dir+(0,))
        self.cell_rads[i] = rad
        self.initCellState(cellState)
        self.set_cells()


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
        kernel_src = resource_string(__name__, 'CLSPP.cl').decode()

        self.program = cl.Program(self.context, kernel_src).build(cache_dir=False)
        # Some kernels that seem like they should be built into pyopencl...
        self.vclearf = ElementwiseKernel(self.context, "float4 *v", "v[i]=0.0", "vecclearf")
        self.vcleari = ElementwiseKernel(self.context, "int *v", "v[i]=0", "veccleari")
        self.vsub = ElementwiseKernel(self.context, "float4 *res, const float4 *in1, const float4 *in2",
                                          "res[i] = in1[i] - in2[i]", "vecsub")
        self.vaddkx = ElementwiseKernel(self.context,
                                            "float4 *res, const float k, const float4 *in1, const float4 *in2",
                                            "res[i] = in1[i] + k*in2[i]", "vecaddkx")
        self.vsubkx = ElementwiseKernel(self.context,
                                            "float4 *res, const float k, const float4 *in1, const float4 *in2",
                                            "res[i] = in1[i] - k*in2[i]", "vecsubkx")

        # cell geometry kernels
        self.calc_cell_area = ElementwiseKernel(self.context, "float* res, float* r, float* l",
                                           "res[i] = 2.f*3.1415927f*r[i]*(2.f*r[i]+l[i])", "cell_area_kern")
        self.calc_cell_vol = ElementwiseKernel(self.context, "float* res, float* r, float* l",
                                          "res[i] = 3.1415927f*r[i]*r[i]*(2.f*r[i]+l[i])", "cell_vol_kern")

        # A dot product as sum of float4 dot products -
        # i.e. like flattening vectors of float4s into big float vectors
        # then computing dot
        # NB. Some openCLs seem not to implement dot(float4,float4) so split
        # into float4's
        #self.vdot = ReductionKernel(self.context, numpy.float32, neutral="0",
        #        reduce_expr="a+b", map_expr="dot(x[i].s0123,y[i].s0123)+dot(x[i].s4567,y[i].s4567)",
        #        arguments="__global float4 *x, __global float4 *y")
    
        # For float4
        self.vdot = ReductionKernel(self.context, numpy.float32, neutral="0",
                reduce_expr="a+b", map_expr="dot(x[i],y[i])",
                arguments="__global float4 *x, __global float4 *y")

        # Add a force to position part of generalised position vector
        self.add_force = ElementwiseKernel(self.context,
                                            "float4 *pos, const float mag, const float4 *dir",
                                            "pos[i].s0123 = pos[i].s0123 + mag*dir[i]", "add_force")

    def init_data(self):
        """Set up the data OpenCL will store on the device."""
        # cell data
        cell_geom = (self.max_cells,)
        self.cell_centers = numpy.zeros(cell_geom, vec.float4)
        self.cell_centers_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_dirs = numpy.zeros(cell_geom, vec.float4)
        self.cell_angs_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_angs = numpy.zeros(cell_geom, numpy.float32)
        self.cell_dirs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.pred_cell_centers = numpy.zeros(cell_geom, vec.float4)
        self.pred_cell_centers_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.pred_cell_dirs = numpy.zeros(cell_geom, vec.float4)
        self.pred_cell_dirs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.cell_rads = numpy.zeros(cell_geom, numpy.float32)
        self.cell_rads_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_sqs = numpy.zeros(cell_geom, numpy.int32)
        self.cell_sqs_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)
        self.cell_n_cts = numpy.zeros(cell_geom, numpy.int32)
        self.cell_n_cts_dev = cl_array.zeros(self.queue, cell_geom, numpy.int32)
        self.cell_dcenters = numpy.zeros(cell_geom, vec.float4)
        self.cell_dcenters_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.avg_neighbour_dir  = numpy.zeros(cell_geom, vec.float4)
        self.avg_neighbour_dir_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)

        self.cell_areas_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32) 
        self.cell_areas_dev[:] = 4 * numpy.pi 
        self.cell_vols_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32) + 4 * numpy.pi / 3
        self.cell_vols_dev[:] = 4 * numpy.pi / 3 
        self.cell_old_vols_dev = self.cell_vols_dev # No growth

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

        self.fr_ents = numpy.zeros(mat_geom, vec.float4)
        self.fr_ents_dev = cl_array.zeros(self.queue, mat_geom, vec.float4)
        self.to_ents = numpy.zeros(mat_geom, vec.float4)
        self.to_ents_dev = cl_array.zeros(self.queue, mat_geom, vec.float4)
        

        # vectors and intermediates
        self.deltap = numpy.zeros(cell_geom, vec.float4)
        self.deltap_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.Mx = numpy.zeros(mat_geom, numpy.float32)
        self.Mx_dev = cl_array.zeros(self.queue, mat_geom, numpy.float32)
        self.BTBx = numpy.zeros(cell_geom, vec.float4)
        self.BTBx_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.Minvx_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)

        # CGS intermediates
        self.p_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.Ap_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.res_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
        self.rhs_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)
    

    def load_from_cellstates(self, cell_states):
        for (cid,cs) in list(cell_states.items()):
            i = cs.idx
            self.cell_centers[i] = tuple(cs.pos)+(0,)
            self.cell_dirs[i] = tuple(cs.dir)+(0,)
            self.cell_rads[i] = cs.radius
        
        self.n_cells = len(cell_states)
        self.set_cells()

    def get_cells(self):
        """Copy cell centers, dirs, lens, and rads from the device."""
        self.cell_centers[0:self.n_cells] = self.cell_centers_dev[0:self.n_cells].get()
        self.cell_dirs[0:self.n_cells] = self.cell_dirs_dev[0:self.n_cells].get()
        self.cell_tos[0:self.n_cells,:] = self.cell_tos_dev[0:self.n_cells,:].get()
        self.avg_neighbour_dir[0:self.n_cells] = self.avg_neighbour_dir_dev[0:self.n_cells].get()
        self.cell_rads[0:self.n_cells] = self.cell_rads_dev[0:self.n_cells].get()
        self.cell_dcenters[0:self.n_cells] = self.cell_dcenters_dev[0:self.n_cells].get()
        self.ct_tos[0:self.n_cts] = self.ct_tos_dev[0:self.n_cts].get()

    def set_cells(self):
        """Copy cell centers, dirs, lens, and rads to the device from local."""
        self.cell_centers_dev[0:self.n_cells].set(self.cell_centers[0:self.n_cells])
        self.cell_dirs_dev[0:self.n_cells].set(self.cell_dirs[0:self.n_cells])
        self.avg_neighbour_dir_dev[0:self.n_cells] = self.avg_neighbour_dir[0:self.n_cells]
        self.cell_rads_dev[0:self.n_cells].set(self.cell_rads[0:self.n_cells])
        self.cell_dcenters_dev[0:self.n_cells].set(self.cell_dcenters[0:self.n_cells])

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
        x_dev = cl_array.zeros(self.queue, (self.n_cells,), vec.float4)
        Ax_dev = cl_array.zeros(self.queue, (self.n_cells,), vec.float4)
        opstring = ''
        for i in range(self.n_cells):
            x = numpy.zeros((self.n_cells,), vec.float4)
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
        if self.printing: print("MTM")
        if self.printing: print(opstring)
        open('CellModeller/Biophysics/BacterialModels/matrix.mat', 'w').write(opstring)


    def dump_cell_data(self, n):
        import pickle
        filename = 'data/data-%04i.pickle'%n
        outfile = open(filename, 'wb')
        data = (self.n_cells,
                self.cell_centers_dev.get(),
                self.cell_dirs_dev.get(),
                self.cell_rads_dev.get(),
                self.parents),
        pickle.dump(data, outfile, protocol=-1)

    def dydt(self):
        self.set_cells()

    def finish(self):
        # pull cells from the device and update simulator
        if self.simulator:
            self.get_cells()
            for id,state in self.simulator.cellStates.items():
                self.updateCellState(state)

    def progress_init(self, dt):
        #self.set_cells()
        # NOTE: by default self.dt=None, and time step == simulator time step (dt) 
        if self.dt:
            self.n_ticks = int(math.ceil(dt/self.dt)) 
        else:
            self.n_ticks = 1 
        #if self.printing: print("n_ticks = %d"%(self.n_ticks))
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
        start = time.time()
        if self.simulator:
            self.get_cells()
            for state in list(self.simulator.cellStates.values()):
                self.updateCellState(state)
        end = time.time()
        if self.printing: print('Updating cell states took ', end-start)

    def step(self, dt):
        """Step forward dt units of time.

        Assumes that:
        cell_centers is up to date when it starts.
        """
        start = time.time()
        if not self.progress_initialised:
            self.progress_init(dt)
        if self.progress():
            self.progress_finalise()
            end = time.time()
            if self.printing: print('step took %g'%(end-start))
            return True
        else:
            return False

    def sub_tick_init(self, dt):
        start = time.time()
        # Compute angle of cell orientation
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

        end = time.time()
        if self.printing: print('sub_tick_init took %g', end-start)

    def tick(self, dt):
        if not self.sub_tick_initialised:
            self.sub_tick_init(dt)
        if self.sub_tick(dt):
            self.sub_tick_finalise(dt)
            return True
        else:
            return False

    def sub_tick(self, dt):
        self.sub_tick_i += 1
        if self.sub_tick_i>self.max_substeps:
            return True
        old_n_cts = self.n_cts
        self.predict()
        # find all contacts
        start = time.time()
        self.find_contacts()
        end1 = time.time()
        if self.printing: print('find_contacts took %g'%(end1-start))
        # place 'backward' contacts in cells
        self.collect_tos()
        end2 = time.time()
        if self.printing: print('collect_tos took %g'%(end2-end1))

        alpha = 10**(self.sub_tick_i)
        new_cts = self.n_cts - old_n_cts
        if (new_cts>0 or self.sub_tick_i==1):
            start = time.time()
            self.build_matrix() # Calculate entries of the matrix
            end = time.time()
            if self.printing: print('build_matrix took %g'%(end-start))
            #print "max cell contacts = %i"%cl_array.max(self.cell_n_cts_dev).get()
            self.CGSSolve(dt, alpha) # invert MTMx to find deltap
            self.add_impulse()
            end = time.time()
            if self.printing: print('CGS solver %g'%(end-start))
            return False
        else:
            return True

    def sub_tick_finalise(self, dt):
        #print "Substeps = %d"%self.sub_tick_i
        self.integrate(dt)
        self.sub_tick_initialised=False

    def initCellState(self, state):
        cid = state.id
        i = state.idx
        state.pos = [self.cell_centers[i][j] for j in range(3)]
        state.dir = [self.cell_dirs[i][j] for j in range(3)]
        state.avg_neighbour_dir = [self.avg_neighbour_dir[i][j] for j in range(3)]
        state.radius = self.cell_rads[i]
        
        pa = numpy.array(state.pos)
        da = numpy.array(state.dir)

    def updateCellState(self, state):
        cid = state.id
        i = state.idx

        state.vel = [self.cell_centers[i][j]-state.pos[j] for j in range(3)]
        state.pos = [self.cell_centers[i][j] for j in range(3)]
        state.dir = [self.cell_dirs[i][j] for j in range(3)]
        state.radius = self.cell_rads[i]
        state.avg_neighbour_dir = [self.avg_neighbour_dir[i][j] for j in range(3)]

        if self.computeNeighbours: #populate cellstate.neighbours
            tos = self.ct_tos[state.idx,:]
            idx2Id = self.simulator.idxToId
            state.neighbours = [idx2Id[to] for to in tos if to>0]
        state.cts = len(state.neighbours)

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
        start = time.time()
        self.sorted_ids.put(numpy.arange(self.n_cells), numpy.argsort(self.cell_sqs[:self.n_cells]))
        self.sorted_ids_dev[0:self.n_cells].set(self.sorted_ids[0:self.n_cells])

        # find the start of each sq in the list of sorted cell ids and send to the device
        sorted_sqs = numpy.sort(self.cell_sqs[:self.n_cells])
        self.sq_inds.put(numpy.arange(self.n_sqs), numpy.searchsorted(sorted_sqs, numpy.arange(self.n_sqs), side='left'))
        self.sq_inds_dev.set(self.sq_inds)
        end = time.time()
        #if self.printing: print('sort_cells took %g'%(end-start))


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
        else:
            centers = self.cell_centers_dev
            dirs = self.cell_dirs_dev

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
                                   numpy.float32(self.Ws),
                                   numpy.float32(self.Wc),
                                   centers.data,
                                   dirs.data,
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
                                   self.ct_overlap_dev.data,
                                   self.avg_neighbour_dir_dev.data).wait()

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
                                  numpy.float32(self.Ws),
                                  numpy.float32(self.Wc),
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

        self.program.calculate_Mx(self.queue,
                                      (self.n_cells,),
                                      None,
                                      numpy.float32(self.gamma_s),
                                      self.cell_dirs_dev.data,
                                      self.cell_rads_dev.data,
                                      x.data,
                                      self.Mx_dev.data).wait()
        # Correct scaling of viscous drag by 1/dt
        self.vaddkx(Ax, 1/dt, Ax, self.Mx_dev).wait()
    

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

        self.add_force(self.rhs_dev, self.Fm, self.cell_dirs_dev).wait()
        #self.add_force(self.rhs_dev, -0.001, self.cell_centers_dev).wait()

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
                if self.printing: print('% 5i'%self.frame_no + '% 6i cells  % 6i cts  % 6i iterations  residual = %f' % (self.n_cells,
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

        if self.frame_no%10==0:
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
                             self.cell_dcenters_dev.data,
                             self.pred_cell_centers_dev.data).wait()

    def integrate(self, dt):
        """Integrates cell centers, dirs, lens for a timestep dt based
        on the current deltap.

        Assumes cell_centers, cell_dirs, cell_lens, cell_rads, and
        deltap are current on the device.

        Calculates new cell_centers, cell_dirs, cell_lens.
        """
        noise = self.rand.normal(self.queue, sigma=np.sqrt(dt), shape=(self.n_cells,), dtype=np.float32)
        if self.simulator.integ:
            signal_gradient = self.simulator.integ.cellSigGradients_dev.data
        else:
            signal_gradient = None
        if self.simulator.sig:
            nSignals = self.simulator.sig.nSignals
        else:
            nSignals = 0
        self.program.integrate(self.queue,
                               (self.n_cells,),
                               None,
                               self.cell_centers_dev.data,
                               self.cell_dirs_dev.data,
                               self.cell_dcenters_dev.data,
                               self.avg_neighbour_dir_dev.data,
                               signal_gradient,
                               noise.data,
                               numpy.float32(self.fcil),
                               numpy.float32(self.ftax),
                               numpy.float32(self.forg),
                               numpy.float32(self.D),
                               numpy.float32(dt),
                               numpy.int32(self.spherical*1),
                               numpy.float32(self.sphere_radius),
                               numpy.int32(nSignals)).wait()

    def add_impulse(self):
        self.program.add_impulse(self.queue, (self.n_cells,), None,
                                 numpy.float32(self.gamma_s),
                                 self.deltap_dev.data,
                                 self.cell_dirs_dev.data,
                                 self.cell_rads_dev.data,
                                 self.cell_dcenters_dev.data).wait()

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
        if self.printing: print("Grid stuff timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
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
        if self.printing: print("Find contacts timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
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
            if self.printing: print("cgs prof: iters=%i, res=%f"%(iters,res))
        t2 = time.clock()
        if self.printing: print("CGS timing for 1000 calls, time per call (s) = %f"%((t2-t1)*0.001))
        open("cgs_prof","a").write( "%i, %i, %i, %f\n"%(self.n_cells,self.n_cts,iters,(t2-t1)*0.001) )


