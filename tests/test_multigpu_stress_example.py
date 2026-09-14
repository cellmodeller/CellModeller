"""Capacity policy checks for the example; no OpenCL installation required."""
import importlib.util
import ast
import io
import os
import random
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    'stress_example_test', Path(__file__).resolve().parents[1] / 'Examples/multigpu_stress.py')
stress = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stress)


class StressExampleTests(unittest.TestCase):
    def test_default_founder_matches_growth_examples(self):
        with patch.dict(os.environ, {}, clear=True):
            cfg = stress.configuration(None)
        self.assertEqual(cfg['initial_cells'], 1)
        self.assertEqual(cfg['growth_rate'], 2.0)
        self.assertEqual(list(stress.founder_geometry(1, random.Random(123))),
                         [((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))])

    def test_gui_model_uses_standard_setup_and_division_budget(self):
        # Execute model callbacks with stand-ins so this remains GPU-independent.
        path = Path(__file__).resolve().parents[1] / 'Examples/multigpu_stress_gui.py'
        tree = ast.parse(path.read_text())
        tree.body = [node for node in tree.body if not isinstance(node, (ast.Import, ast.ImportFrom))]
        calls = []
        class Physics:
            def __init__(self, sim, **kw):
                calls.append(kw)
        namespace = dict(random=random, capacity_estimate=lambda *args: 123456,
                         host_available_bytes=lambda: 2**30,
                         CLBacterium=Physics,
                         ModuleRegulator=lambda sim: 'regulator')
        exec(compile(tree, str(path), 'exec'), namespace)
        sim = SimpleNamespace(is_gui=False, stepNum=0, CLWorkStats={},
                              init=lambda *args: calls.append(args),
                              addCell=lambda **kw: calls.append(kw),
                              setSaveOutput=lambda enabled: None)
        namespace['setup'](sim)
        self.assertIsInstance(calls[1][0], Physics)
        self.assertEqual(calls[1][1:], ('regulator', None, None))
        self.assertEqual(calls[0]['max_cells'], 123456)
        self.assertEqual(calls[0]['max_sqs'], 4 * 123456)
        self.assertEqual(calls[-1], dict(cellType=0, pos=(0, 0, 0)))
        namespace['_capacity'] = 3
        cells = {i: SimpleNamespace(volume=4, targetVol=3) for i in range(2)}
        namespace['update'](cells)
        self.assertEqual(sum(c.divideFlag for c in cells.values()), 1)
        self.assertTrue(all(c.growthRate == 0 for c in cells.values()))

    def test_simultaneous_divisions_cannot_exceed_capacity(self):
        stress._cfg = dict(max_cells=10, growth_rate=1.0, report_every=10)
        stress._sim = SimpleNamespace(is_gui=False)
        cells = {i: SimpleNamespace(volume=5, targetVol=4, divideFlag=True) for i in range(8)}
        stress.update(cells)
        self.assertEqual(sum(cell.divideFlag for cell in cells.values()), 2)
        self.assertTrue(all(cell.growthRate == 0 for cell in cells.values()))
        cells.update({i: SimpleNamespace(volume=5, targetVol=4, divideFlag=True) for i in (8, 9)})
        stress.update(cells)
        self.assertFalse(any(cell.divideFlag for cell in cells.values()))

    def test_growth_continues_when_capacity_remains(self):
        stress._cfg = dict(max_cells=10, growth_rate=2.0, report_every=10)
        stress._sim = SimpleNamespace(is_gui=False)
        cell = SimpleNamespace(volume=2, targetVol=4, divideFlag=True)
        stress.update({0: cell})
        self.assertFalse(cell.divideFlag)
        self.assertEqual(cell.growthRate, 2.0)

    def test_auto_capacity_respects_host_and_execution_mode(self):
        device = SimpleNamespace(global_mem_size=2**30, max_mem_alloc_size=2**28)
        sim = SimpleNamespace(CLDevices=[device, device], clDeviceWeights=[0.5, 0.5],
                              CLDevicePool=SimpleNamespace(partitioned=True))
        partitioned = stress.capacity_estimate(sim, 32, 4, 0.5, 2**40)
        sim.CLDevicePool.partitioned = False
        replicated = stress.capacity_estimate(sim, 32, 4, 0.5, 2**40)
        self.assertGreater(partitioned, replicated)
        self.assertLessEqual(stress.capacity_estimate(sim, 32, 4, 0.5, 2**20),
                             int(2**19 / (8192 + 768 * 32 + 32 * 4)))

    def test_configuration_rejects_overflow_and_small_capacity(self):
        for env in ({'CM_STRESS_MAX_CELLS': '0'},
                    {'CM_STRESS_MAX_CELLS': str(2**30)},
                    {'CM_STRESS_MAX_CONTACTS': '4'},
                    {'CM_STRESS_GROWTH_RATE': 'nan'}):
            with self.subTest(env=env), patch.dict(os.environ, env, clear=True):
                with self.assertRaises(ValueError):
                    stress.configuration(None)

    def test_report_records_interval_counts_without_mutating_simulator(self):
        stress._cfg = dict(max_cells=100)
        stress._previous = {'physics.find_contacts': [10, 10]}
        stress._last_report = 0
        sim = SimpleNamespace(stepNum=2, cellStates={1: None},
                              CLWorkStats={'physics.find_contacts': [20, 21]},
                              CLMemoryStats={'physics.find_contacts': [{'bytes': 32}, {'bytes': 40}]})
        with patch('builtins.print'), patch.object(stress, 'host_available_bytes', return_value=1000):
            result = stress.report(sim, io.StringIO())
        self.assertEqual(result['work_items']['physics.find_contacts'], [10, 11])
        self.assertEqual(sim.CLWorkStats['physics.find_contacts'], [20, 21])


if __name__ == '__main__':
    unittest.main()
