"""Dependency-free scheduling tests: python -m unittest discover -s tests -v."""
import ast
import importlib.util
from pathlib import Path
import re
import sys
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class Buffer:
    def __init__(self, context, flags, size):
        self.size = size
        self.data = bytearray(size)


class Queue:
    def __init__(self, label):
        self.context = 'context'
        self.label = label
        self.pending = []

    def finish(self):
        while self.pending:
            self.pending.pop(0)()

    def flush(self):
        pass


class Event:
    def __init__(self, queue):
        self.queue = queue

    def wait(self):
        self.queue.finish()


def copy(queue, dest, src, src_offset=0, dst_offset=0, byte_count=None):
    size = src.size if byte_count is None else byte_count
    def action():
        if src_offset+size > src.size or dst_offset+size > dest.size:
            raise AssertionError('Out-of-bounds transfer')
        dest.data[dst_offset:dst_offset+size] = src.data[src_offset:src_offset+size]
    queue.pending.append(action)
    return Event(queue)


fake_cl = types.SimpleNamespace(Buffer=Buffer, enqueue_copy=copy,
                                enqueue_marker=Event,
                                mem_flags=types.SimpleNamespace(READ_WRITE=1),
                                device_type=types.SimpleNamespace(GPU=4))
spec = importlib.util.spec_from_file_location('multi_gpu_under_test', ROOT/'CellModeller/MultiGPU.py')
mg = importlib.util.module_from_spec(spec)
with patch.dict(sys.modules, {'pyopencl': fake_cl}):
    spec.loader.exec_module(mg)


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.devices = [types.SimpleNamespace(type=t) for t in (2, 4, 4)]

    def test_legacy_cpu_and_gpu(self):
        self.assertEqual(mg.device_indices(self.devices, None), [0])
        self.assertEqual(mg.device_indices(self.devices, None, 2), [2])

    def test_auto_and_single_device_fallback(self):
        self.assertEqual(mg.device_indices(self.devices, 'auto'), [1, 2])
        self.assertEqual(mg.device_indices(self.devices[:1], 'auto'), [0])
        self.assertEqual(mg.device_indices(self.devices[1:2], 'auto'), [0])

    def test_explicit_order(self):
        self.assertEqual(mg.device_indices(self.devices, [2, 1]), [2, 1])

    def test_invalid_selections(self):
        for selection in ([], [-1], [3], [True], [1.0], [1, 1], [0, 1], 'all', 1):
            with self.subTest(selection=selection), self.assertRaises(ValueError):
                mg.device_indices(self.devices, selection)

    def test_weights(self):
        self.assertEqual(mg.normalized_weights(None, 2), [.5, .5])
        self.assertEqual(mg.normalized_weights([3, 1], 2), [.75, .25])
        self.assertEqual(mg.normalized_weights([1e308, 1e308], 2), [.5, .5])

    def test_invalid_weights(self):
        for weights in ([], [1], [1, 0], [-1, 2], [float('nan'), 1], [float('inf'), 1]):
            with self.subTest(weights=weights), self.assertRaises(ValueError):
                mg.normalized_weights(weights, 2)

    def test_partition_coverage(self):
        for count in range(100):
            for weights in ([1, 1], [1, 3, 2], [1000, 1], [1, 1000]):
                parts = list(mg.partitions(count, mg.normalized_weights(weights, len(weights))))
                cells = [i for start, size in parts for i in range(start, start+size)]
                self.assertEqual(cells, list(range(count)))


class DispatchTests(unittest.TestCase):
    def fixture(self, name, count=11, weights=(1, 2, 1), threshold=1):
        queues = [Queue(str(i)) for i in range(len(weights))]
        first, stride_arg, overlap = mg.ContactProgram.layouts[name]
        widths = [4] + [3*b for b in (4, 4, 4, 16, 16, 4, 4)] + ([12] if overlap else [])
        args = [0]*first
        args[stride_arg] = 3
        # Full geometry and contact state, including inactive rows.
        args[first-1] = Buffer('context', 1, 64)
        args[first-1].data[:] = bytes(range(64))
        for width in widths:
            buf = Buffer('context', 1, (count+2)*width)
            buf.data[:] = bytes([7])*buf.size
            args.append(buf)
        calls = []

        def kernel(queue, global_size, local_size, *kernel_args, **kwargs):
            start = kwargs.get('global_offset', (0,))[0]
            size = global_size[0]
            calls.append((queue.label, start, size))
            def action():
                # Read remote geometry (outside the worker's owned range), and
                # update existing outputs. This catches missing input snapshots.
                self.assertEqual(kernel_args[first-1].data[-1], 63)
                for index, width in enumerate(widths, first):
                    for cell in range(start, start+size):
                        offset = cell*width
                        self.assertEqual(kernel_args[index].data[offset], 7)
                        kernel_args[index].data[offset:offset+width] = bytes([cell+20])*width
            queue.pending.append(action)
            return Event(queue)

        program = types.SimpleNamespace(**{name: kernel, 'other': kernel})
        wrapped = mg.ContactProgram(program, queues,
                                    mg.normalized_weights(weights, len(weights)), threshold)
        return queues, args, widths, calls, wrapped, first

    def test_all_contact_outputs_merge_and_preserve_inactive_rows(self):
        for name in mg.ContactProgram.layouts:
            with self.subTest(kernel=name):
                queues, args, widths, calls, wrapped, first = self.fixture(name)
                getattr(wrapped, name)(queues[0], (11,), None, *args).wait()
                self.assertEqual(calls, [('0', 0, 3), ('1', 3, 5), ('2', 8, 3)])
                for index, width in enumerate(widths, first):
                    expected = b''.join(bytes([i+20])*width for i in range(11)) + bytes([7])*width*2
                    self.assertEqual(args[index].data, expected)

    def test_small_populations_and_single_device(self):
        for weights, threshold in (((1, 1), 12), ((1,), 1)):
            queues, args, _, calls, wrapped, _ = self.fixture('find_contacts', weights=weights, threshold=threshold)
            wrapped.find_contacts(queues[0], (11,), None, *args).wait()
            self.assertEqual(calls, [('0', 0, 11)])

    def test_more_devices_than_cells(self):
        queues, args, _, calls, wrapped, _ = self.fixture('find_contacts', count=1)
        wrapped.find_contacts(queues[0], (1,), None, *args).wait()
        self.assertEqual(calls, [('1', 0, 1)])

    def test_unaudited_kernel_passthrough(self):
        _, _, _, _, wrapped, _ = self.fixture('find_contacts')
        self.assertIs(wrapped.other, wrapped.program.other)

    def test_reject_unsupported_launch(self):
        queues, args, _, _, wrapped, _ = self.fixture('find_contacts')
        with self.assertRaises(ValueError):
            wrapped.find_contacts(queues[0], (11,), (1,), *args)

    def test_output_layout_matches_actual_kernel_signatures(self):
        source = (ROOT/'CellModeller/Biophysics/BacterialModels/CLBacterium.cl').read_text()
        for name, (first, stride_arg, overlap) in mg.ContactProgram.layouts.items():
            signature = re.search(r'__kernel void '+name+r'\((.*?)\)\s*\{', source, re.S).group(1)
            arguments = signature.split(',')
            self.assertTrue(arguments[stride_arg].strip().endswith('max_contacts'))
            expected = ['n_cts', 'frs', 'tos', 'dists', 'pts', 'norms', 'reldists', 'stiff']
            if overlap:
                expected.append('overlap')
            self.assertEqual([a.split('*')[-1].strip() for a in arguments[first:]], expected)


class SimulatorConfigurationTests(unittest.TestCase):
    def setUp(self):
        # Exercise the actual initialization method without importing the
        # unrelated numerical dependencies required by Simulator.py.
        tree = ast.parse((ROOT/'CellModeller/Simulator.py').read_text())
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef))
        method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'init_cl')
        self.devices = [types.SimpleNamespace(type=4, name='GPU %d' % i) for i in range(2)]
        platform = types.SimpleNamespace(name='test', get_devices=lambda: self.devices)
        self.contexts = []
        def context(**kwargs):
            self.contexts.append(kwargs)
            return kwargs
        fake = types.SimpleNamespace(get_platforms=lambda: [platform], Context=context,
                                     CommandQueue=lambda ctx, device: (ctx, device),
                                     context_properties=types.SimpleNamespace(PLATFORM=1))
        namespace = {'cl': fake, '__package__': 'CellModeller'}
        exec(compile(ast.Module(body=[method], type_ignores=[]), 'Simulator.py', 'exec'), namespace)
        self.initialize = namespace['init_cl']
        self.sim = types.SimpleNamespace()
        self.modules = patch.dict(sys.modules, {'CellModeller.MultiGPU': mg})
        self.modules.start()
        self.addCleanup(self.modules.stop)

    def test_primary_and_worker_queues_match_requested_order(self):
        with patch('builtins.print'):
            self.initialize(self.sim, 0, 0, [1, 0], [3, 1], 12)
        self.assertEqual(self.sim.CLDevices, self.devices[::-1])
        self.assertEqual(self.sim.CLQueue, self.sim.CLQueues[0])
        self.assertEqual(self.sim.clDeviceWeights, [.75, .25])
        self.assertEqual(self.sim.clMultiGPUMinCells, 12)

    def test_legacy_initialization(self):
        with patch('builtins.print'):
            self.initialize(self.sim, 0, 1)
        self.assertEqual(self.sim.CLDevices, [self.devices[1]])
        self.assertEqual(len(self.sim.CLQueues), 1)

    def test_bad_platform_and_threshold_fail_before_context_creation(self):
        for platform, threshold in ((-1, 1), (1, 1), (True, 1), (0, 0), (0, 1.5), (0, True)):
            with self.subTest(platform=platform, threshold=threshold), self.assertRaises(ValueError):
                self.initialize(self.sim, platform, 0, min_cells=threshold)
        self.assertEqual(self.contexts, [])


if __name__ == '__main__':
    unittest.main()
