"""Resident CG execution against an independent dense NumPy reference.

The fake runtime executes the resident kernel operations and records transfers;
real OpenCL coverage lives in test_multi_gpu_hardware.py.
"""
import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

try:
    import numpy as np
except ImportError:
    np = None


class Event:
    def wait(self):
        pass


class Array:
    def __init__(self, value):
        self.value = value
        self.data = self
    def __getitem__(self, key):
        return Array(self.value[key])
    def get(self):
        return self.value.copy()


def copy(queue, target, source, byte_count=None, device_offset=0, **kwargs):
    left = target.value if isinstance(target, Array) else target
    right = source.value if isinstance(source, Array) else source
    a, b = left.view(np.uint8).reshape(-1), right.view(np.uint8).reshape(-1)
    count = len(b) if byte_count is None else byte_count
    a[device_offset:device_offset+count] = b[:count]
    return Event()


class Program:
    fail = False
    retrievals = 0
    def __getattribute__(self, name):
        if name.startswith('rd_') or name == 'calculate_Mx':
            type(self).retrievals += 1
        return object.__getattribute__(self, name)
    def __init__(self, *args):
        pass
    def build(self):
        return self
    def rd_transpose(self, q, shape, local, offsets, contacts, signs, f, t, bx, result):
        for i in range(shape[0]):
            result.value[i] = 0
            for k in range(offsets.value[i], offsets.value[i+1]):
                c = contacts.value[k]
                result.value[i] += (f.value[c] if signs.value[k] > 0 else -t.value[c])*bx.value[c]
    def rd_bx(self, q, shape, local, fr, to, valid, f, t, p, bx):
        if self.fail:
            raise RuntimeError('injected device failure')
        for i in range(shape[0]):
            bx.value[i] = (np.dot(f.value[i], p.value[fr.value[i]]) -
                           (np.dot(t.value[i], p.value[to.value[i]]) if to.value[i] >= 0 else 0)) if valid.value[i] else 0
    def calculate_Mx(self, q, shape, local, mu, gamma, dirs, lens, rads, p, mx):
        for i in range(shape[0]):
            # Fixture uses the x axis: independent diagonal inertia expression.
            l = lens.value[i]+2*rads.value[i]
            mx.value[i] = p.value[i]*np.array([mu*l]*3+[0, mu*l**3/12, mu*l**3/12, gamma, 0])
    def rd_regularize(self, q, shape, local, scale, ap, mx):
        ap.value += scale*mx.value
    def rd_update(self, q, shape, local, alpha, x, r, p, ap):
        x.value += alpha*p.value[:shape[0]]
        r.value -= alpha*ap.value
    def rd_direction(self, q, shape, local, beta, p, r):
        p.value[:shape[0]] = r.value+beta*p.value[:shape[0]]
    def rd_gather(self, q, shape, local, indices, p, packed):
        packed.value[:] = p.value[indices.value]


@unittest.skipIf(np is None, 'Requires NumPy')
class ResidentTests(unittest.TestCase):
    def setUp(self):
        path = Path(__file__).resolve().parents[1] / 'CellModeller/ResidentSolver.py'
        spec = importlib.util.spec_from_file_location('resident_under_test', path)
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        arrays = types.ModuleType('pyopencl.array')
        arrays.to_device = lambda q, v: Array(v.copy())
        arrays.zeros = lambda q, n, dtype: Array(np.zeros(n, dtype))
        arrays.empty = lambda q, n, dtype: Array(np.empty(n, dtype))
        arrays.vec = types.SimpleNamespace(float8=np.dtype((np.float32, 8)))
        cl = types.ModuleType('pyopencl')
        cl.array, cl.Program, cl.enqueue_copy = arrays, Program, copy
        reduction = types.ModuleType('pyopencl.reduction')
        reduction.ReductionKernel = lambda *a, **kw: lambda x, y, **kw: Array(np.asarray(np.sum(x.value*y.value), np.float32))
        context = patch.dict(sys.modules, {'pyopencl': cl, 'pyopencl.array': arrays,
                                           'pyopencl.reduction': reduction})
        context.start()
        self.addCleanup(context.stop)
        Program.fail = False

    def fixture(self, n=6):
        stride = 8
        owned = [np.arange(n//2, dtype=np.int32), np.arange(n//2, n, dtype=np.int32)]
        queue = types.SimpleNamespace(context=0, finish=lambda: None, flush=lambda: None)
        pool = types.SimpleNamespace(queues=[queue, queue], ownership=lambda count: owned,
                                     clear_cache=lambda: None, check_memory=lambda *args: None,
                                     record=lambda *args: None)
        phys = types.SimpleNamespace(n_cells=n, max_contacts=stride, gamma=2., muA=1., cgs_tol=1e-5)
        def field(name, value):
            setattr(phys, name, types.SimpleNamespace(array=np.asarray(value)))
        fr = np.repeat(np.arange(n, dtype=np.int32), stride)
        to = np.repeat((np.arange(n, dtype=np.int32)+1)%n, stride)
        entries = np.tile(np.eye(8, dtype=np.float32), (n, 1))
        rhs = np.random.default_rng(1).normal(size=n*stride).astype(np.float32)
        field('cell_n_cts_dev', np.full(n, stride, np.int32))
        field('n_cell_tos_dev', np.full(n, stride, np.int32))
        field('cell_tos_dev', np.concatenate([np.arange(((i-1)%n)*stride, ((i-1)%n+1)*stride, dtype=np.int32) for i in range(n)]))
        for name, value in [('ct_frs_dev', fr), ('ct_tos_dev', to), ('fr_ents_dev', entries),
                            ('to_ents_dev', .1*entries), ('ct_reldists_dev', rhs),
                            ('cell_dirs_dev', np.tile(np.array([1,0,0,0], np.float32), (n,1))),
                            ('cell_lens_dev', np.full(n, 2., np.float32)),
                            ('cell_rads_dev', np.full(n, .5, np.float32)),
                            ('deltap_dev', np.full((n,8), -123., np.float32))]:
            field(name, value)
        b = np.zeros((n*8, n*8))
        for cell in range(n):
            for k in range(8):
                b[cell*8+k, cell*8+k] = 1
                b[cell*8+k, ((cell+1)%n)*8+k] = -.1
        m = np.tile([3,3,3,0,27/12,27/12,2,0], n)
        expected = np.linalg.solve(b.T@b+np.diag(m/2), b.T@rhs).reshape(n,8)
        return pool, phys, expected

    def test_resident_iterations_match_dense_reference(self):
        pool, phys, expected = self.fixture()
        solver = self.module.ResidentCG(pool)
        retrievals = Program.retrievals
        iterations, residual = solver.solve(phys)
        self.assertEqual(Program.retrievals, retrievals)
        self.assertGreater(iterations, 0)
        self.assertLess(residual, phys.cgs_tol)
        np.testing.assert_allclose(phys.deltap_dev.array, expected, rtol=1e-4, atol=1e-5)
        self.assertEqual(solver.stats['solution_bytes'], phys.n_cells*32)
        # Ring partitions exchange four boundary vectors per iteration, each way.
        self.assertEqual(solver.stats['halo_bytes'], iterations*4*32*2)
        self.assertTrue(solver.stats['converged'])

    def test_failure_does_not_commit_solution(self):
        pool, phys, _ = self.fixture()
        solver = self.module.ResidentCG(pool)
        Program.fail = True
        with self.assertRaisesRegex(RuntimeError, 'injected'):
            solver.solve(phys)
        self.assertTrue(np.all(phys.deltap_dev.array == -123))

    def test_zero_rhs_requires_no_halo_exchange(self):
        pool, phys, _ = self.fixture()
        phys.ct_reldists_dev.array[:] = 0
        solver = self.module.ResidentCG(pool)
        self.assertEqual(solver.solve(phys), (0, 0.0))
        self.assertEqual(solver.stats['halo_bytes'], 0)
        self.assertTrue(np.all(phys.deltap_dev.array == 0))

    def test_preflight_failure_prevents_solution_commit(self):
        pool, phys, _ = self.fixture()
        def fail(*args):
            raise MemoryError('resident capacity')
        pool.check_memory = fail
        with self.assertRaises(MemoryError):
            self.module.ResidentCG(pool).solve(phys)
        self.assertTrue(np.all(phys.deltap_dev.array == -123))


if __name__ == '__main__':
    unittest.main()
