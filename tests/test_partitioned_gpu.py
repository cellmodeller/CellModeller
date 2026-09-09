"""Partition ownership and capacity tests; NumPy required, no OpenCL hardware."""
import importlib.util
from pathlib import Path
import re
import sys
import types
import unittest
from unittest.mock import patch

try:
    import numpy as np
except ImportError:
    np = None

ROOT = Path(__file__).resolve().parents[1]


class FakeBuffer:
    allocations = []
    def __init__(self, context, flags, size=None, hostbuf=None):
        self.context = context
        self.data = (np.ascontiguousarray(hostbuf).view(np.uint8).reshape(-1).copy()
                     if hostbuf is not None else np.zeros(size, np.uint8))
        self.size = self.data.nbytes
        self.allocations.append((context, self.size))


class Event:
    def wait(self):
        pass


class Queue:
    def __init__(self, device, memory=2**30):
        self.context = device
        self.device = types.SimpleNamespace(global_mem_size=memory, max_mem_alloc_size=memory//2, name=str(device))
    def finish(self):
        pass
    def flush(self):
        pass


def copy(queue, target, source):
    if source.context != queue.context:
        raise AssertionError('Cross-context copy')
    target.view(np.uint8).reshape(-1)[:] = source.data
    return Event()


class Program:
    sources = []
    fail_device = None
    def __init__(self, context, source):
        self.context = context
        self.sources.append(source)
    def build(self, **kwargs):
        return self
    def advance(self, queue, shape, local, inp, out, owners, error, in_ids, in_count, in_width, out_ids, out_count, out_width):
        for buf in (inp, out, owners, error, in_ids, out_ids):
            if buf.context != queue.context:
                raise AssertionError('Buffer assigned to the wrong device')
        if queue.context == self.fail_device:
            error.data.view(np.uint32)[0] = 1
        source = inp.data.view(np.float32)
        target = out.data.view(np.float32)
        source_keys, target_keys = in_ids.data.view(np.int32), out_ids.data.view(np.int32)
        for gid in owners.data.view(np.int32):
            target[np.searchsorted(target_keys, gid)] = 2*source[np.searchsorted(source_keys, gid)]
        return Event()


fake_cl = types.SimpleNamespace(Buffer=FakeBuffer, Program=Program, enqueue_copy=copy,
                                mem_flags=types.SimpleNamespace(READ_WRITE=1, COPY_HOST_PTR=2))


@unittest.skipIf(np is None, 'Partitioned tests require NumPy')
class PartitionedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load under the package so that relative imports use the real scheduler.
        spec = importlib.util.spec_from_file_location('CellModeller.partitioned_under_test', ROOT/'CellModeller/PartitionedGPU.py')
        cls.module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {'pyopencl': fake_cl, spec.name: cls.module}):
            spec.loader.exec_module(cls.module)

    def setUp(self):
        cl_patch = patch.dict(sys.modules, {'pyopencl': fake_cl})
        cl_patch.start()
        self.addCleanup(cl_patch.stop)
        self.p = self.module
        FakeBuffer.allocations = []
        Program.sources = []
        Program.fail_device = None

    def test_arrays_allocate_no_device_memory_and_preserve_slice_copy(self):
        array = self.p.PartitionedArrays.zeros(Queue(0), (100000, 24), np.float32)
        self.assertEqual(FakeBuffer.allocations, [])
        array[10:12].set(np.ones((2, 24), np.float32))
        self.assertEqual(array.get()[10:12].sum(), 48)
        target = self.p.PartitionedArrays.zeros(Queue(0), (2, 24), np.float32)
        self.p.copy_array(Queue(0), target.data, array[10:12].data).wait()
        np.testing.assert_array_equal(target.get(), np.ones((2, 24)))
        self.assertEqual(FakeBuffer.allocations, [])

    def test_spatial_partition_remaps_interleaved_global_ids(self):
        xyz = np.array([[0, 0, 0], [100, 0, 0], [1, 0, 0], [101, 0, 0]], np.float32)
        order = self.p.spatial_order(xyz, [.5, .5])
        np.testing.assert_array_equal(order, [0, 2, 1, 3])
        self.assertEqual(sorted(order.tolist()), list(range(4)))

    def test_empty_and_weighted_spatial_partitions(self):
        for count in range(30):
            xyz = np.zeros((count, 3), np.float32)
            order = self.p.spatial_order(xyz, [.25, .5, .25])
            self.assertEqual(sorted(order.tolist()), list(range(count)))

    def fixture(self, devices=2, count=100):
        pool = self.p.PartitionedPool([Queue(i) for i in range(devices)], [1/devices]*devices, 10000)
        source = '__kernel void advance(__global const float *inp, __global float *out) { int i=get_global_id(0); out[i]=2*inp[i]; }'
        program = self.p.PartitionedProgram(source, pool, {'advance': lambda a: [(1, 4, 0)]}, family='elementwise')
        inp = self.p.PartitionedArray(np.arange(count, dtype=np.float32))
        out = self.p.PartitionedArray(np.zeros(count, np.float32))
        return pool, program, inp, out

    def test_no_primary_replica_and_device_footprint_shrinks(self):
        two, program, inp, out = self.fixture(2, 1000)
        program.advance(two.queues[0], (1000,), None, inp.data, out.data).wait()
        np.testing.assert_array_equal(out.get(), inp.get()*2)
        usage_two = [s['bytes'] for s in two.memory_stats['advance']]
        self.assertTrue(all(size < inp.nbytes for _, size in FakeBuffer.allocations))
        self.assertEqual(two.stats['advance'], [500, 500])
        FakeBuffer.allocations = []
        four, program, inp, out = self.fixture(4, 1000)
        program.advance(four.queues[0], (1000,), None, inp.data, out.data).wait()
        self.assertTrue(all(s['bytes'] < usage_two[0] for s in four.memory_stats['advance']))
        self.assertEqual(four.stats['advance'], [250]*4)
        self.assertEqual({context for context, _ in FakeBuffer.allocations}, {0, 1, 2, 3})
        self.assertTrue(all(size <= 1000 for _, size in FakeBuffer.allocations))

    def test_errors_do_not_commit_partial_worker_results(self):
        pool, program, inp, out = self.fixture()
        Program.fail_device = 1
        with self.assertRaisesRegex(RuntimeError, 'missing halo'):
            program.advance(pool.queues[0], (100,), None, inp.data, out.data)
        np.testing.assert_array_equal(out.get(), np.zeros(100))

    def test_capacity_is_checked_before_any_allocation(self):
        pool, program, inp, out = self.fixture()
        for queue in pool.queues:
            queue.device.global_mem_size = 10
        with self.assertRaises(MemoryError):
            program.advance(pool.queues[0], (100,), None, inp.data, out.data)
        self.assertEqual(FakeBuffer.allocations, [])

    def test_owned_output_commit_preserves_unowned_rows(self):
        source = self.p.HostBuffer(np.arange(40, dtype=np.float32))
        region = self.p.Region([1, 3, 7], 4)
        packed = region.pack(source, 4)
        packed[:] = 0
        region.commit(source, 4, packed, [3])
        expected = np.arange(40, dtype=np.float32)
        expected[12:16] = 0
        np.testing.assert_array_equal(source.array, expected)

    def test_kernel_rewriting_covers_every_pointer_access_and_preserves_helpers(self):
        paths = [ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl']
        paths += list((ROOT/'CellModeller/Integration').glob('*.cl'))
        for path in paths:
            source = path.read_text()
            transformed = self.p.partition_source(source)
            parts = self.p.kernel_parts(source)
            self.assertEqual(transformed.count('__kernel void '), len(parts), path)
            self.assertEqual(transformed.count('cm_owned[get_global_id(0)]'), len(parts), path)
            for name, (_, _, _, body, args) in parts.items():
                for arg in args:
                    if arg.pointer:
                        self.assertIn('cm_'+arg.name+'_ids', transformed)
            self.assertIn('atomic_or(error, 1u)', transformed)
        nested = self.p.rewrite_accesses('out[indices[i]] = source[indices[i]];', ['out', 'source', 'indices'])
        self.assertEqual(nested.count('cm_lookup('), 4)
        expression = self.p.rewrite_accesses('rates = specRate+specbase;', ['specRate'])
        self.assertIn('specRate + cm_lookup(', expression)

    def test_contact_halo_is_local_and_grid_index_is_sparse(self):
        from CellModeller.MultiGPU import PHYSICS_LAYOUTS
        source = (ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        specs = self.p.kernel_parts(source)['find_contacts'][4]
        count, stride = 100, 4
        scalars = dict(max_cells=count, n_cells=count, grid_x_min=0, grid_x_max=count,
                       grid_y_min=0, grid_y_max=1, n_sqs=count, max_contacts=stride)
        args = []
        for arg in specs:
            if not arg.pointer:
                args.append(scalars[arg.name])
            else:
                dtype = np.int32 if arg.dtype == 'int' else np.float32
                data = np.zeros(count*stride*4, dtype)
                if arg.name in ('sqs', 'sorted_ids', 'sq_inds'):
                    data[:count] = np.arange(count)
                args.append(self.p.HostBuffer(data))
        regions = self.p.plan_regions('find_contacts', specs, args, np.arange(20, 30),
                                      PHYSICS_LAYOUTS['find_contacts'](args), 'physics')
        regions = {arg.name: region for arg, region in zip(specs, regions) if arg.pointer}
        np.testing.assert_array_equal(regions['centers'].ids, np.arange(19, 31))
        np.testing.assert_array_equal(regions['frs'].ids, np.arange(20, 30))
        np.testing.assert_array_equal(regions['sorted_ids'].ids, np.arange(19, 31))
        self.assertLess(len(regions['sq_inds'].ids), count)

    def test_transpose_exchange_contains_only_referenced_remote_contacts(self):
        from CellModeller.MultiGPU import PHYSICS_LAYOUTS
        source = (ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        specs = self.p.kernel_parts(source)['calculate_BTBx'][4]
        stride = 8
        counts = np.zeros(100, np.int32); counts[10:12] = 1
        tos = np.zeros(800, np.int32); tos[80], tos[88] = 3, 723
        values = {'max_contacts': stride, 'n_cts': counts, 'n_cell_tos': counts, 'cell_tos': tos,
                  'fr_ents': np.zeros(800*8, np.float32), 'to_ents': np.zeros(800*8, np.float32),
                  'Bx': np.zeros(800, np.float32), 'BTBx': np.zeros(800, np.float32)}
        args = [self.p.HostBuffer(values[a.name]) if a.pointer else values[a.name] for a in specs]
        regions = self.p.plan_regions('calculate_BTBx', specs, args, np.array([10, 11]),
                                      PHYSICS_LAYOUTS['calculate_BTBx'](args), 'physics')
        result = {a.name: r for a, r in zip(specs, regions) if a.pointer}
        np.testing.assert_array_equal(result['to_ents'].ids, [3, 723])
        np.testing.assert_array_equal(result['Bx'].ids, [3, 80, 88, 723])
        self.assertEqual(result['to_ents'].width, 1)

    def test_signal_grid_halo_excludes_unused_grid_nodes(self):
        from CellModeller.MultiGPU import SIGNAL_LAYOUTS
        source = (ROOT/'CellModeller/Integration/CLCrankNicIntegrator.cl').read_text()
        specs = self.p.kernel_parts(source)['setCellSignals'][4]
        values = dict(numSignals=2, gridTotalSize=1000, gridDimx=10, gridDimy=10, gridDimz=10,
                      indices=np.tile(np.arange(8, dtype=np.int32), 100), weights=np.zeros(800, np.float32),
                      grid=np.zeros(2000, np.float32), levels=np.zeros(200, np.float32))
        args = [self.p.HostBuffer(values[a.name]) if a.pointer else values[a.name] for a in specs]
        regions = self.p.plan_regions('setCellSignals', specs, args, np.arange(10), SIGNAL_LAYOUTS['setCellSignals'](args), 'signal')
        result = {a.name: r for a, r in zip(specs, regions) if a.pointer}
        np.testing.assert_array_equal(result['grid'].ids, list(range(8))+list(range(1000, 1008)))
        self.assertEqual(len(result['levels'].ids), 10)
        self.assertEqual(result['levels'].width, 2)

    def test_physics_allocator_never_allocates_full_state_on_a_gpu(self):
        import ast
        tree = ast.parse((ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.py').read_text())
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef))
        method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'init_data')
        vector = lambda n: np.dtype([(str(i), np.float32) for i in range(n)])
        namespace = {'numpy': np, 'vec': types.SimpleNamespace(float4=vector(4), float8=vector(8))}
        exec(compile(ast.Module(body=[method], type_ignores=[]), 'CLBacterium.py', 'exec'), namespace)
        model = types.SimpleNamespace(max_cells=1000, max_contacts=24, max_planes=1,
                                      max_spheres=1, max_sqs=200, device_arrays=self.p.PartitionedArrays, queue=Queue(0))
        namespace['init_data'](model)
        self.assertEqual(FakeBuffer.allocations, [])
        for name, value in vars(model).items():
            if name.endswith('_dev'):
                self.assertIsInstance(value, self.p.PartitionedArray, name)
        self.assertEqual(model.fr_ents_dev.size, 24000)

    def test_partitioned_reduction_only_uploads_local_views(self):
        pool = self.p.PartitionedPool([Queue(0), Queue(1)], [.5, .5], 10000)
        pool.cell_order = np.array([0, 2, 4, 6, 1, 3, 5, 7], np.int32)
        values = self.p.PartitionedArray(np.arange(12, dtype=np.float32)[2:10])
        uploads, calls = [], []
        def upload(queue, data):
            uploads.append((queue.context, data.copy()))
            return types.SimpleNamespace(context=queue.context, value=data)
        def factory(context):
            def reduce(x, y, queue):
                self.assertEqual(x.context, context)
                self.assertEqual(queue.context, context)
                calls.append(context)
                return types.SimpleNamespace(get=lambda: np.asarray(np.dot(x.value, y.value), np.float32))
            return reduce
        arrays = types.SimpleNamespace(to_device=upload)
        with patch.object(fake_cl, 'array', arrays, create=True), patch.dict(sys.modules, {'pyopencl.array': arrays}):
            reduction = self.p.PartitionedReduction(factory(0), pool, np.float32, 'dot', factory)
            result = reduction(values, values).get()
        self.assertEqual(calls, [0, 1])
        self.assertTrue(all(len(data) == 4 for _, data in uploads))
        self.assertEqual(result.item(), np.dot(values.array, values.array))
        self.assertEqual(pool.stats['dot'], [4, 4])

    def test_full_integer_region_limits_are_validated_before_casting(self):
        with self.assertRaises(ValueError):
            self.p.Region([2**32])
        with self.assertRaises(ValueError):
            self.p.Region([-1])
        with self.assertRaises(ValueError):
            self.p.Region([0], 0)

    def test_contact_capacity_errors_are_guarded_in_partitioned_source(self):
        source = (ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        result = self.p.partition_source(source)
        self.assertGreaterEqual(result.count('k >= max_contacts'), 4)
        self.assertIn('atomic_or(cm_error, 2u)', result)

    def test_reverse_contact_scan_uploads_only_occupied_target_slots(self):
        from CellModeller.MultiGPU import PHYSICS_LAYOUTS
        source = (ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        specs = self.p.kernel_parts(source)['collect_tos'][4]
        count, stride = 100, 24
        scalars = dict(max_cells=count, n_cells=count, grid_x_min=0, grid_x_max=count,
                       grid_y_min=0, grid_y_max=1, n_sqs=count, max_contacts=stride)
        args = []
        for arg in specs:
            if not arg.pointer:
                args.append(scalars[arg.name])
            else:
                values = np.zeros(count*stride, np.int32)
                if arg.name in ('sqs', 'sorted_ids', 'sq_inds'):
                    values[:count] = np.arange(count)
                elif arg.name == 'n_cts':
                    values[:count] = 1
                args.append(self.p.HostBuffer(values))
        regions = self.p.plan_regions('collect_tos', specs, args, np.arange(20, 30),
                                      PHYSICS_LAYOUTS['collect_tos'](args), 'physics')
        by_name = {arg.name: region for arg, region in zip(specs, regions) if arg.pointer}
        np.testing.assert_array_equal(by_name['tos'].ids, np.arange(19, 31)*stride)
        self.assertEqual(by_name['tos'].width, 1)
        self.assertEqual(len(by_name['frs'].ids), 0)


if __name__ == '__main__':
    unittest.main()
