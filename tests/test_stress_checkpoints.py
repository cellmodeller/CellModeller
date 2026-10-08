"""Checkpoint durability and runner continuation without OpenCL hardware."""
import importlib.util
import json
from pathlib import Path
import random
import signal
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('stress_checkpoint_tests',
    Path(__file__).resolve().parents[1] / 'Examples/multigpu_stress.py')
stress = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stress)


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.np = SimpleNamespace(random=SimpleNamespace(
            get_state=lambda: ('numpy-state',), set_state=lambda state: None))
        self.model = SimpleNamespace(_cfg={'max_cells': 10}, _rng=random.Random(42))
        self.sim = SimpleNamespace(dt=.025, stepNum=20, cellStates={1: 'cell'},
            lineage={2: 1}, _next_id=4, _next_idx=1,
            integ=SimpleNamespace(levels=[.2, .4]), CLWorkStats={})

    def save(self, path):
        with patch.dict(sys.modules, {'numpy': self.np}):
            stress.save_checkpoint(path, self.sim, self.model, 15, 100)

    def test_roundtrip_restores_rng_counters_and_species(self):
        p = self.root/'latest.pickle'
        self.save(p)
        expected = [self.model._rng.random() for _ in range(5)]
        data = stress.read_checkpoint(p)
        def load(data):
            self.sim.cellStates = data['cellStates']
            self.sim.stepNum = data['stepNum']
            self.sim.integ.levels = data['specData']
        self.sim.loadFromPickle = load
        self.sim._next_id = 999
        with patch.dict(sys.modules, {'numpy': self.np}):
            stress.restore_checkpoint(self.sim, self.model, data)
        self.assertEqual([self.model._rng.random() for _ in range(5)], expected)
        self.assertEqual(self.sim._next_id, 4)
        self.assertEqual(self.sim.integ.nCells, 1)
        self.assertEqual(self.sim.integ.levels, [.2, .4])
        self.assertEqual(data['reached'], 15)
        self.assertEqual(data['config']['max_cells'], 10)

    def test_failed_write_preserves_previous_checkpoint(self):
        p = self.root/'latest.pickle'
        self.save(p)
        original = p.read_bytes()
        with patch.object(stress.pickle, 'dump', side_effect=OSError('disk full')):
            with self.assertRaises(OSError): self.save(p)
        self.assertEqual(p.read_bytes(), original)
        self.assertEqual(list(self.root.glob('*.tmp')), [])

    def test_source_mismatch_rejected(self):
        p = self.root/'latest.pickle'
        self.save(p)
        with patch.object(stress, 'checkpoint_signature', return_value='changed'):
            with self.assertRaisesRegex(ValueError, 'source differs'):
                stress.read_checkpoint(p)

    def test_runner_safe_stop_and_resume_preserve_hold_progress(self):
        handlers = {}
        def install(sig, fn):
            old = handlers.get(sig, signal.SIG_DFL)
            handlers[sig] = fn
            return old
        stop_on_step = [2]
        class FakeSimulator:
            def __init__(obj, path, dt, **kw):
                obj.dt = dt; obj.stepNum = 0
                cfg = stress.configuration(None)
                obj.cellStates = {1: 'cell'}; obj.lineage = {}
                obj._next_id = 2; obj._next_idx = 1
                obj.integ = SimpleNamespace(levels=[.2])
                obj.CLDevices = [SimpleNamespace(name='fake', global_mem_size=100, driver_version='fake')]
                obj.clDeviceWeights = [1.]; obj.CLWorkStats = {}
                obj.module = SimpleNamespace(_cfg=cfg, _rng=random.Random(cfg['seed']),
                    validate_state=lambda sim: None, report=lambda sim, stream: None)
            def step(obj):
                obj.stepNum += 1
                if obj.stepNum == stop_on_step[0]: handlers[signal.SIGINT](signal.SIGINT, None)
            def loadFromPickle(obj, data):
                obj.cellStates = data['cellStates']; obj.stepNum = data['stepNum']
                obj.integ.levels = data['specData']
        p = self.root/'latest.pickle'
        modules = {'numpy': self.np, 'CellModeller.Simulator': SimpleNamespace(Simulator=FakeSimulator)}
        with patch.dict(sys.modules, modules), patch.object(stress.signal, 'signal', side_effect=install):
            with patch.object(sys, 'argv', ['stress', '--max-cells', '1', '--hold-steps', '4',
                        '--checkpoint', str(p), '--log', str(self.root/'first.jsonl')]):
                stress.main()
            d = stress.read_checkpoint(p)
            self.assertEqual(d['stepNum'], 2)
            self.assertEqual(d['reached'], 1)
            self.assertEqual(json.loads((self.root/'first.jsonl').read_text().splitlines()[-1])['status'], 'stopped_checkpointed')
            stop_on_step[0] = -1
            with patch.object(sys, 'argv', ['stress', '--resume', str(p), '--log', str(self.root/'second.jsonl')]):
                stress.main()
            d = stress.read_checkpoint(p)
            self.assertEqual(d['stepNum'], 5)  # completes remaining hold, not a new full hold
            self.assertEqual(d['config']['max_cells'], 1)  # not CLI default 100000
            self.assertEqual(json.loads((self.root/'second.jsonl').read_text().splitlines()[-1])['status'], 'target_held')


if __name__ == '__main__':
    unittest.main()
