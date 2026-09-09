"""Scheduling and ownership checks for all distributed simulation stages."""
import ast
from pathlib import Path
import re
import struct
import sys
import types
import unittest
from unittest.mock import patch

from test_multi_gpu import mg, Buffer, Queue, Event, fake_cl, ROOT


class StageTests(unittest.TestCase):
    def test_every_physics_and_integration_kernel_has_a_layout(self):
        cases = [
            ('Biophysics/BacterialModels/CLBacterium.cl', mg.PHYSICS_LAYOUTS),
            ('Integration/CLEulerIntegrator.cl', mg.EULER_LAYOUTS),
            ('Integration/CLCrankNicIntegrator.cl', mg.SIGNAL_LAYOUTS),
            ('Integration/CLEulerSigIntegrator.cl', mg.SIGNAL_LAYOUTS),
        ]
        for path, layouts in cases:
            source = (ROOT/'CellModeller'/path).read_text()
            names = set(re.findall(r'__kernel void (\w+)\(', source))
            self.assertEqual(names, set(layouts), path)
            for name in names:
                signature = re.search(r'__kernel void '+name+r'\((.*?)\)\s*\{', source, re.S).group(1)
                declarations = signature.split(',')
                outputs = layouts[name]([3]*len(declarations))
                for index, width, offset in outputs:
                    self.assertIn('*', declarations[index], (name, index))
                    self.assertGreater(width, 0)
                    self.assertEqual(offset, 0)

    def test_two_dimensional_kernels_split_cells_and_keep_contact_dimension(self):
        queues = [Queue(str(i)) for i in range(3)]
        pool = mg.DevicePool(queues, mg.normalized_weights(None, 3), 1)
        output = Buffer('context', 1, 8*5*4)
        calls = []
        def kernel(queue, shape, local, arg, global_offset):
            calls.append((queue.label, shape, global_offset))
            start = global_offset[0]
            def action():
                for cell in range(start, start+shape[0]):
                    for contact in range(shape[1]):
                        struct.pack_into('i', arg.data, (cell*5+contact)*4, cell*100+contact)
            queue.pending.append(action)
            return Event(queue)
        program = mg.DistributedProgram(types.SimpleNamespace(matrix=kernel), pool,
                                         {'matrix': mg.rows((0, 20))})
        program.matrix(queues[0], (8, 5), None, output).wait()
        self.assertEqual(calls, [('0', (3, 5), (0, 0)), ('1', (3, 5), (3, 0)), ('2', (2, 5), (6, 0))])
        self.assertEqual(pool.stats['matrix'], [15, 15, 10])
        self.assertEqual(list(struct.unpack('40i', output.data)),
                         [cell*100+contact for cell in range(8) for contact in range(5)])

    def test_equal_devices_receive_equal_work_for_divisible_populations(self):
        for devices in range(2, 9):
            weights = mg.normalized_weights(None, devices)
            self.assertEqual([size for _, size in mg.partitions(100*devices, weights)],
                             [100]*devices)

    def test_replica_reused_refreshed_and_aliases_preserved(self):
        queues = [Queue('0'), Queue('1')]
        pool = mg.DevicePool(queues, [.5, .5], 1)
        buffer = Buffer('context', 1, 16)
        identities = []
        def kernel(queue, shape, local, out, inp, global_offset):
            self.assertIs(out, inp)
            if queue is queues[1]:
                identities.append(id(out))
            def action():
                for cell in range(global_offset[0], global_offset[0]+shape[0]):
                    value = struct.unpack_from('i', inp.data, cell*4)[0]
                    struct.pack_into('i', out.data, cell*4, value+1)
            queue.pending.append(action)
            return Event(queue)
        program = mg.DistributedProgram(types.SimpleNamespace(update=kernel), pool, {'update': mg.rows((0, 4))})
        for value in (10, 50):
            buffer.data[:] = struct.pack('4i', *([value]*4))
            program.update(queues[0], (4,), None, buffer, buffer).wait()
            self.assertEqual(struct.unpack('4i', buffer.data), (value+1,)*4)
        self.assertEqual(identities[0], identities[1])

    def test_cache_budget_and_explicit_release(self):
        pool = mg.DevicePool([Queue('0'), Queue('1')], [.5, .5], 1, cache_bytes=32)
        a, b, large = [Buffer('context', 1, n) for n in (24, 24, 40)]
        first = pool.replica(1, a)
        self.assertIs(first, pool.replica(1, a))
        pool.replica(1, b)
        self.assertEqual(pool.cache_sizes, [0, 24])
        self.assertIsNot(pool.replica(1, a), first)
        pool.replica(1, large)
        self.assertEqual(pool.cache_sizes, [0, 24])
        pool.clear_cache()
        self.assertEqual(pool.cache_sizes, [0, 0])
        self.assertFalse(pool.caches[1])

    def test_offset_output_merges_only_the_view(self):
        queues = [Queue('0'), Queue('1')]
        pool = mg.DevicePool(queues, [.5, .5], 1)
        out = Buffer('context', 1, 32)
        out.data[:] = bytes([9])*32
        def kernel(queue, shape, local, arg, global_offset):
            def action():
                start = 8+global_offset[0]*4
                arg.data[start:start+shape[0]*4] = bytes([1])*shape[0]*4
            queue.pending.append(action)
            return Event(queue)
        program = mg.DistributedProgram(types.SimpleNamespace(view=kernel), pool,
                                         {'view': lambda args: [(0, 4, 8)]})
        program.view(queues[0], (4,), None, out).wait()
        self.assertEqual(out.data, bytes([9])*8+bytes([1])*16+bytes([9])*8)

    def test_no_launch_for_empty_work_and_invalid_outputs_fail_early(self):
        queue = Queue('0')
        calls = []
        kernel = lambda *a, **k: calls.append(a)
        pool = mg.DevicePool([queue, Queue('1')], [.5, .5], 1)
        program = mg.DistributedProgram(types.SimpleNamespace(k=kernel), pool, {'k': mg.rows((0, 4))})
        buffer = Buffer('context', 1, 4)
        program.k(queue, (0,), None, buffer).wait()
        with self.assertRaises(ValueError):
            program.k(queue, (2,), None, buffer)
        self.assertEqual(calls, [])


class FakeDtype:
    kind = 'f'
    itemsize = 4


class Array:
    def __init__(self, queue, shape, dtype, data=None, offset=0):
        self.queue = queue
        self.size = shape[0]
        self.dtype = dtype
        self.offset = offset
        self.base_data = data
        self.nbytes = self.size*dtype.itemsize
        self.flags = types.SimpleNamespace(c_contiguous=True)

    def values(self):
        return struct.unpack_from('%df' % self.size, self.base_data.data, self.offset)


class Scalar:
    def __init__(self, value):
        self.value = value

    def item(self):
        return self.value


class ReductionTests(unittest.TestCase):
    def test_reductions_launch_all_partitions_before_readback_and_sum_globally(self):
        queues = [Queue(str(i)) for i in range(3)]
        pool = mg.DevicePool(queues, mg.normalized_weights(None, 3), 1)
        buffer = Buffer('context', 1, 10*4)
        buffer.data[:] = struct.pack('10f', *range(10))
        view = Array(queues[0], (6,), FakeDtype(), data=buffer, offset=8)
        calls = []
        def reduction(x, y, queue=None):
            calls.append((queue.label, x.size))
            self.assertEqual(x.offset, y.offset)
            class Partial:
                def get(inner_self):
                    self.assertEqual(len(calls), 3)
                    return Scalar(sum(a*b for a, b in zip(x.values(), y.values())))
            return Partial()
        numpy = types.SimpleNamespace(dtype=lambda d: d, asarray=lambda value, dtype: Scalar(value))
        arrays = types.SimpleNamespace(Array=Array, to_device=lambda q, value: types.SimpleNamespace(get=lambda: value))
        with patch.object(fake_cl, 'array', arrays, create=True), patch.dict(sys.modules, {'numpy': numpy, 'pyopencl': fake_cl, 'pyopencl.array': arrays}):
            wrapped = mg.DistributedReduction(reduction, pool, FakeDtype(), 'dot')
            result = wrapped(view, view).get().item()
        self.assertEqual(result, sum(i*i for i in range(2, 8)))
        self.assertEqual(calls, [('0', 2), ('1', 2), ('2', 2)])
        self.assertEqual(pool.stats['dot'], [2, 2, 2])

    def test_below_threshold_uses_original_reduction(self):
        pool = mg.DevicePool([Queue('0'), Queue('1')], [.5, .5], 10)
        array = types.SimpleNamespace(size=2)
        original = lambda *args, **kwargs: 'original result'
        arrays = types.SimpleNamespace()
        with patch.object(fake_cl, 'array', arrays, create=True), patch.dict(sys.modules, {'numpy': types.SimpleNamespace(), 'pyopencl': fake_cl,
                                      'pyopencl.array': arrays}):
            wrapped = mg.DistributedReduction(original, pool, FakeDtype(), 'dot')
            self.assertEqual(wrapped(array, array), 'original result')
        self.assertEqual(pool.stats['dot'], [2, 0])


class ElementwiseTests(unittest.TestCase):
    def test_generated_offsets_and_aliased_vector_update(self):
        queues = [Queue('0'), Queue('1')]
        pool = mg.DevicePool(queues, [.5, .5], 1)
        buffer = Buffer('context', 1, 8*4)
        buffer.data[:] = struct.pack('8f', *range(8))
        view = Array(queues[0], (4,), FakeDtype(), data=buffer, offset=8)
        sources = []
        class Program:
            def __init__(inner_self, context, source):
                sources.append(source)
            def build(inner_self, **kwargs):
                return inner_self
            def add(inner_self, queue, shape, local, output, out_offset, scale, inp, in_offset, global_offset):
                def action():
                    for i in range(global_offset[0], global_offset[0]+shape[0]):
                        value = struct.unpack_from('f', inp.data, (in_offset+i)*4)[0]
                        old = struct.unpack_from('f', output.data, (out_offset+i)*4)[0]
                        struct.pack_into('f', output.data, (out_offset+i)*4, old+scale*value)
                queue.pending.append(action)
                return Event(queue)
        numpy = types.SimpleNamespace(uint64=int, float32=float, int32=int)
        with patch.object(fake_cl, 'Program', Program, create=True), patch.dict(sys.modules, {'numpy': numpy}):
            wrapped = mg.DistributedElementwise(None, pool, 'float *out, const float k, const float *inp',
                                                 'out[i] += k*inp[i]', 'add')
            wrapped(view, 2, view).wait()
        self.assertIn('out_data + out_offset', sources[0])
        self.assertEqual(struct.unpack('8f', buffer.data), (0, 1, 6, 9, 12, 15, 6, 7))
        self.assertEqual(pool.stats['add'], [2, 2])


if __name__ == '__main__':
    unittest.main()
