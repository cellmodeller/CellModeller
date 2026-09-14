"""Optional real-kernel comparison; skips when two GPUs are unavailable."""
from pathlib import Path
import re
import unittest

try:
    import numpy as np
    import pyopencl as cl
except ImportError:
    np = cl = None


@unittest.skipIf(cl is None or np is None, 'Requires NumPy and PyOpenCL')
class HardwareSimulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            platforms = cl.get_platforms()
        except cl.Error as error:
            raise unittest.SkipTest(str(error))
        for platform in platforms:
            devices = [d for d in platform.get_devices() if d.type & cl.device_type.GPU]
            if len(devices) >= 2:
                cls.platform_index = platforms.index(platform)
                available = platform.get_devices()
                cls.device_indices = [available.index(d) for d in devices[:2]]
                cls.context = cl.Context(devices=devices[:2])
                cls.queues = [cl.CommandQueue(cls.context, device=d) for d in devices[:2]]
                break
        else:
            raise unittest.SkipTest('Requires two GPUs on one OpenCL platform')
        from CellModeller.MultiGPU import ContactProgram
        cls.wrapper = ContactProgram
        cls.source = (Path(__file__).resolve().parents[1]/
                      'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        cls.program = cl.Program(cls.context, cls.source).build()

    def arguments(self, name):
        # All cells share a grid square. Adjacent overlapping cells straddle
        # the split, exercising contacts that need geometry from another GPU.
        count, capacity, contacts = 7, 9, 32
        scalars = dict(max_cells=capacity, n_cells=count, max_contacts=contacts,
                       n_planes=1, n_spheres=1, n_sqs=1,
                       grid_x_min=0, grid_x_max=1, grid_y_min=0, grid_y_max=1)
        signature = re.search(r'__kernel void '+name+r'\((.*?)\)\s*\{', self.source, re.S).group(1)
        args, arrays = [], {}
        for declaration in signature.split(','):
            arg_name = declaration.split()[-1].lstrip('*')
            if '*' not in declaration:
                args.append(np.int32(scalars[arg_name]))
                continue
            dtype = np.int32 if 'int*' in declaration else np.float32
            vector = 4 if 'float4*' in declaration else 1
            rows = capacity*contacts if arg_name in ('frs', 'tos', 'dists', 'pts', 'norms', 'reldists', 'stiff', 'overlap') else capacity
            array = np.zeros((rows, vector), dtype=dtype)
            if arg_name == 'centers':
                array[:, 0] = np.arange(rows)*0.7
                array[:, 2] = 0.25
            elif arg_name == 'dirs':
                array[:, 1] = 1
            elif arg_name in ('lens', 'plane_coeffs', 'sphere_coeffs', 'sphere_norms'):
                array[:] = 1
            elif arg_name == 'rads':
                array[:] = 0.5
            elif arg_name == 'sorted_ids':
                array[:, 0] = np.arange(rows)
            elif arg_name == 'plane_norms':
                array[:, 2] = 1
            elif arg_name == 'sphere_pts':
                array[:, 2] = -1
            elif arg_name == 'sphere_rads':
                array[:] = 1
            buffer = cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=array)
            args.append(buffer)
            arrays[len(args)-1] = array
        return count, args, arrays

    def test_contacts_match_single_gpu_including_repeated_updates(self):
        for name, (first, _, _) in self.wrapper.layouts.items():
            with self.subTest(kernel=name):
                count, single, arrays = self.arguments(name)
                _, multi, _ = self.arguments(name)
                wrapped = self.wrapper(self.program, self.queues, [0.4, 0.6], 1)
                # Second run updates existing contacts, testing replica freshness.
                for repeat in range(2):
                    getattr(self.program, name)(self.queues[0], (count,), None, *single).wait()
                    getattr(wrapped, name)(self.queues[0], (count,), None, *multi).wait()
                    for index in range(first, len(single)):
                        expected = np.empty_like(arrays[index])
                        actual = np.empty_like(expected)
                        cl.enqueue_copy(self.queues[0], expected, single[index]).wait()
                        cl.enqueue_copy(self.queues[0], actual, multi[index]).wait()
                        if expected.dtype.kind == 'i':
                            np.testing.assert_array_equal(actual, expected)
                        else:
                            np.testing.assert_allclose(actual, expected, rtol=1e-4, atol=1e-5)

    def test_all_physics_kernels_against_single_device(self):
        from CellModeller.MultiGPU import DevicePool, DistributedProgram, PHYSICS_LAYOUTS
        from CellModeller.PartitionedGPU import PartitionedPool, PartitionedProgram, HostBuffer
        partition_pool = PartitionedPool(self.queues, [.5, .5], 1)
        partitioned = PartitionedProgram(self.source, partition_pool, PHYSICS_LAYOUTS, 'physics.')
        count, capacity, contacts = 8, 10, 4
        scalars = dict(max_cells=capacity, n_cells=count, max_contacts=contacts,
                       n_planes=0, n_spheres=0, n_sqs=1, grid_x_min=0,
                       grid_x_max=1, grid_y_min=0, grid_y_max=1,
                       grid_spacing=100, muA=1, gamma=10)
        for name, layout in PHYSICS_LAYOUTS.items():
            with self.subTest(kernel=name):
                signature = re.search(r'__kernel void '+name+r'\((.*?)\)\s*\{', self.source, re.S).group(1)
                single, multi, arrays = [], [], {}
                for index, declaration in enumerate(signature.split(',')):
                    arg_name = declaration.split()[-1].lstrip('*')
                    if '*' not in declaration:
                        value = np.float32(scalars[arg_name]) if 'float' in declaration else np.int32(scalars[arg_name])
                        single.append(value)
                        multi.append(value)
                        continue
                    dtype = np.int32 if 'int*' in declaration else np.float32
                    vector = 8 if 'float8*' in declaration else 4 if 'float4*' in declaration else 1
                    array = np.zeros((capacity*contacts, vector), dtype=dtype)
                    if arg_name in ('n_cts', 'n_cell_tos'):
                        array[:count] = 1
                    elif arg_name in ('frs', 'tos', 'cell_tos'):
                        for cell in range(count):
                            value = cell if arg_name == 'frs' else (cell+1)%count if arg_name == 'tos' else ((cell-1)%count)*contacts
                            array[cell*contacts, 0] = value
                        if arg_name == 'tos' and name == 'calculate_Bx':
                            array[0] = -1  # Regression: boundary destination must not be read.
                    elif arg_name == 'sorted_ids':
                        array[:count, 0] = np.arange(count)
                    elif dtype == np.float32:
                        array[:] = 0.1
                        if arg_name in ('dirs', 'pred_dirs'):
                            array[:] = 0
                            array[:, 1] = 1
                        elif arg_name in ('lens', 'rads', 'stiff'):
                            array[:] = 1
                        elif arg_name in ('centers', 'pred_centers'):
                            array[:] = 0
                            array[:, 0] = np.arange(array.shape[0])*3
                    arrays[index] = array
                    for args in (single, multi):
                        args.append(cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=array))
                pool = DevicePool(self.queues, [.5, .5], 1)
                wrapped = DistributedProgram(self.program, pool, PHYSICS_LAYOUTS)
                shape = (count, contacts) if name in ('build_matrix', 'calculate_Bx') else (count,)
                getattr(self.program, name)(self.queues[0], shape, None, *single).wait()
                getattr(wrapped, name)(self.queues[0], shape, None, *multi).wait()
                host_args = [HostBuffer(arrays[i].copy()) if i in arrays else value for i, value in enumerate(single)]
                getattr(partitioned, name)(self.queues[0], shape, None, *host_args).wait()
                for index, _, _ in layout(single):
                    expected, actual = np.empty_like(arrays[index]), np.empty_like(arrays[index])
                    cl.enqueue_copy(self.queues[0], expected, single[index]).wait()
                    cl.enqueue_copy(self.queues[0], actual, multi[index]).wait()
                    if expected.dtype.kind == 'i':
                        np.testing.assert_array_equal(actual, expected)
                    else:
                        np.testing.assert_allclose(actual, expected, rtol=1e-4, atol=1e-5)
                    if expected.dtype.kind == 'i':
                        np.testing.assert_array_equal(host_args[index].array, expected)
                    else:
                        np.testing.assert_allclose(host_args[index].array, expected, rtol=1e-4, atol=1e-5)
                self.assertEqual(pool.stats[name][0], pool.stats[name][1])

    def test_elementwise_and_reduction_views(self):
        import pyopencl.array as array
        from pyopencl.elementwise import ElementwiseKernel
        from pyopencl.reduction import ReductionKernel
        from CellModeller.MultiGPU import DevicePool, DistributedElementwise, DistributedReduction
        pool = DevicePool(self.queues, [.5, .5], 1)
        source = np.arange(20, dtype=np.float32)
        values = array.to_device(self.queues[0], source)
        arguments = 'float *res, const float k, const float *x'
        operation = 'res[i] += k*x[i]'
        original = ElementwiseKernel(self.context, arguments, operation, 'update_values')
        elementwise = DistributedElementwise(original, pool, arguments, operation, 'update_values')
        elementwise(values[2:18], np.float32(2), values[2:18]).wait()
        source[2:18] *= 3
        np.testing.assert_allclose(values.get(), source)
        reduction = ReductionKernel(self.context, np.float32, neutral='0', reduce_expr='a+b',
                                    map_expr='x[i]*y[i]', arguments='__global float *x, __global float *y')
        distributed = DistributedReduction(reduction, pool, np.float32, 'dot')
        actual = distributed(values[2:18], values[2:18]).get()
        np.testing.assert_allclose(actual, np.dot(source[2:18], source[2:18]), rtol=1e-5)
        self.assertEqual(pool.stats['dot'], [8, 8])

    def test_species_and_signal_kernels(self):
        from CellModeller.MultiGPU import DevicePool, DistributedProgram, EULER_LAYOUTS, SIGNAL_LAYOUTS
        root = Path(__file__).resolve().parents[1]/'CellModeller/Integration'
        for file, layouts in [('CLEulerIntegrator.cl', EULER_LAYOUTS),
                              ('CLEulerSigIntegrator.cl', SIGNAL_LAYOUTS),
                              ('CLCrankNicIntegrator.cl', SIGNAL_LAYOUTS)]:
            source = (root/file).read_text()
            if file == 'CLEulerIntegrator.cl':
                source = source % 'rates[0] = species[0]*0.25f; rates[1] = species[1]*0.5f;'
            else:
                source = source % dict(nSignals=1, sigKernel='rates[0] = species[0]*0.1f;',
                                       specKernel='rates[0] = species[0]*0.25f; rates[1] = species[1]*0.5f; species[0] += 0.01f;')
            program = cl.Program(self.context, source).build()
            from CellModeller.PartitionedGPU import PartitionedPool, PartitionedProgram, HostBuffer
            partition_pool = PartitionedPool(self.queues, [.5, .5], 1)
            family = 'euler' if file == 'CLEulerIntegrator.cl' else 'signal'
            partitioned = PartitionedProgram(source, partition_pool, layouts, family=family)
            for name, layout in layouts.items():
                with self.subTest(file=file, kernel=name):
                    signature = re.search(r'__kernel void '+name+r'\((.*?)\)\s*\{', source, re.S).group(1)
                    single, multi, arrays = [], [], {}
                    for index, declaration in enumerate(signature.split(',')):
                        arg_name = declaration.split()[-1].lstrip('*')
                        if '*' not in declaration:
                            value = {'numSpecies': 2, 'numSignals': 1, 'gridTotalSize': 64,
                                     'gridDimx': 4, 'gridDimy': 4, 'gridDimz': 4,
                                     'gridOrigx': 0, 'gridOrigy': 0, 'gridOrigz': 0,
                                     'gridSizex': 1, 'gridSizey': 1, 'gridSizez': 1,
                                     'gridVolume': 1, 'dt': 0.1}[arg_name]
                            value = np.float32(value) if 'float' in declaration else np.int32(value)
                            single.append(value)
                            multi.append(value)
                            continue
                        dtype = np.int32 if 'int*' in declaration else np.float32
                        array = np.ones((128, 4 if 'float4*' in declaration else 1), dtype=dtype)
                        if dtype == np.int32:
                            array[:] = 0
                        elif arg_name == 'weights':
                            array[:] = 0.125
                        arrays[index] = array
                        for args in (single, multi):
                            args.append(cl.Buffer(self.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=array))
                    pool = DevicePool(self.queues, [.5, .5], 1)
                    wrapped = DistributedProgram(program, pool, layouts)
                    getattr(program, name)(self.queues[0], (8,), None, *single).wait()
                    getattr(wrapped, name)(self.queues[0], (8,), None, *multi).wait()
                    host_args = [HostBuffer(arrays[i].copy()) if i in arrays else value for i, value in enumerate(single)]
                    getattr(partitioned, name)(self.queues[0], (8,), None, *host_args).wait()
                    for index, _, _ in layout(single):
                        expected, actual = np.empty_like(arrays[index]), np.empty_like(arrays[index])
                        cl.enqueue_copy(self.queues[0], expected, single[index]).wait()
                        cl.enqueue_copy(self.queues[0], actual, multi[index]).wait()
                        np.testing.assert_allclose(actual, expected, rtol=1e-4, atol=1e-5)
                        np.testing.assert_allclose(host_args[index].array, expected, rtol=1e-4, atol=1e-5)
                    self.assertEqual(pool.stats[name], [4, 4])

    def test_complete_simulation_solver_and_species_with_division(self):
        import random
        from CellModeller.Simulator import Simulator
        model = r"""
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator

def setup(sim):
    physics = CLBacterium(sim, max_cells=32, max_contacts=16, max_sqs=256,
                          jitter_z=False, printing=False, cgs_tol=1e-5)
    regulator = ModuleRegulator(sim)
    integrator = CLEulerIntegrator(sim, 1, 32)
    sim.init(physics, regulator, None, integrator)
    for i in range(6):
        sim.addCell(pos=(i*0.95, 0, 0), dir=(0, 1, 0), length=1.5)
    physics.addPlane((0, 0, -0.49), (0, 0, 1), 0.5)

def init(cell):
    cell.growthRate = 0.1
    cell.species[:] = 1.0 + 0.1*cell.id

def update(cells):
    pass

def specRateCL():
    return 'rates[0] = 0.1f * species[0];'
"""
        single = Simulator('multi_gpu_validation_single', 0.01, moduleStr=model,
                           clPlatformNum=self.platform_index, clDeviceNum=self.device_indices[0])
        multi = Simulator('multi_gpu_validation_multi', 0.01, moduleStr=model,
                          clPlatformNum=self.platform_index, clDeviceNums=self.device_indices,
                          clMultiGPUMinCells=1)
        from CellModeller.PartitionedGPU import PartitionedArray
        self.assertIsInstance(multi.phys.cell_centers_dev, PartitionedArray)
        self.assertIsInstance(multi.integ.specLevel_dev, PartitionedArray)
        self.assertNotEqual(multi.CLQueues[0].context, multi.CLQueues[1].context)
        self.assertIsNotNone(multi.phys.resident_solver)
        for step in range(4):
            if step == 2:
                for sim in (single, multi):
                    random.seed(123)
                    np.random.seed(123)
                    sim.divide(sim.cellStates[min(sim.cellStates)])
            single.step()
            multi.step()
            self.assertEqual(sorted(single.cellStates), sorted(multi.cellStates))
            for cid in single.cellStates:
                left, right = single.cellStates[cid], multi.cellStates[cid]
                np.testing.assert_allclose(right.pos, left.pos, rtol=2e-3, atol=2e-4)
                np.testing.assert_allclose(right.dir, left.dir, rtol=2e-3, atol=2e-4)
                np.testing.assert_allclose(right.length, left.length, rtol=2e-3, atol=2e-4)
                np.testing.assert_allclose(right.species, left.species, rtol=2e-3, atol=2e-4)
        resident = multi.CLResidentSolverStats
        self.assertIn('residual', resident)
        self.assertEqual(resident['solution_bytes'], len(multi.cellStates)*32)
        self.assertEqual(sum(resident['owned_cells']), len(multi.cellStates))
        self.assertTrue(np.isfinite(resident['residual']))
        for stage in ('physics.calculate_Bx', 'physics.calculate_BTBx', 'physics.dot',
                      'vecaddkx', 'CLEulerIntegrator.speciesRates', 'physics.integrate'):
            self.assertIn(stage, multi.CLWorkStats)
            self.assertTrue(all(count > 0 for count in multi.CLWorkStats[stage]), stage)


if __name__ == '__main__':
    unittest.main()
