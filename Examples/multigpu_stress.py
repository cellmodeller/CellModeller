"""Growing-colony capacity test: mechanics, CG solver, division and species.

GUI: load this file; select all desired GPUs, equal weights and partitioned mode.
Default: grow one founder toward 100,000 cells, then stop growth/division.
For larger GUI runs set CM_STRESS_MAX_CELLS=auto before launching the GUI.

Headless (recommended for capacity measurements, from the repository root):
    python Examples/multigpu_stress.py --devices auto --max-cells auto
    python Examples/multigpu_stress.py --devices 0,1 --max-cells 200000
    python Examples/multigpu_stress.py --devices 0 --max-cells 200000
    python Examples/multigpu_stress.py --devices 0,1,2,3 --checkpoint latest.pickle
    python Examples/multigpu_stress.py --devices 0,1,2,3 --resume latest.pickle --log resumed.jsonl

With --checkpoint, Ctrl+C/SIGTERM requests a stop after the current step and an
atomic save. Wait for CHECKPOINT saved and stopped_checkpointed before ending
the session. Periodic checkpoints default to every 100 completed steps. Resume
restores the workload and RNG states; device selection remains configurable.
Only load trusted checkpoints. See Examples/multigpu_checkpointing.md.

Run --help for controls. 'auto' is a conservative memory-based target estimate,
NOT a measured maximum or a promise of full GPU utilization. Increase explicit
targets across separate runs to find the practical limit. Host staging, Python,
contact halos and transfers can bottleneck before GPU memory or compute fills.
No dummy allocations or extra arithmetic are used to inflate utilization.

GUI configuration uses CM_STRESS_ plus MAX_CELLS, INITIAL_CELLS, MAX_CONTACTS,
MAX_SQS, SPECIES, SEED, GROWTH_RATE, REPORT_EVERY, MEMORY_FRACTION (see --help).
This example covers mechanics and intracellular species, not signal diffusion.
Run tests/test_multi_gpu_hardware.py for numerical comparison coverage.
"""

import argparse
import hashlib
import pickle
import signal
import tempfile
import json
import math
import os
from pathlib import Path
import random
import sys
import time


def host_available_bytes():
    """Linux available RAM, bounded by container memory headroom when exposed."""
    limits = []
    try:
        fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
        limits.append(int(fields['MemAvailable'].split()[0]) * 1024)
    except (OSError, KeyError, ValueError):
        pass
    for root, limit_name, usage_name in (
        ('/sys/fs/cgroup', 'memory.max', 'memory.current'),
        ('/sys/fs/cgroup/memory', 'memory.limit_in_bytes', 'memory.usage_in_bytes'),
    ):
        try:
            limit = int((Path(root) / limit_name).read_text())
            used = int((Path(root) / usage_name).read_text())
            limits.append(max(0, limit - used))
        except (OSError, ValueError):
            pass
    return min(limits) if limits else None


def capacity_estimate(sim, contacts, species, fraction, available):
    """Planning heuristic including host copies, cell objects, and halo reserve."""
    if available is None:
        raise ValueError('Cannot estimate available host RAM; set an explicit --max-cells')
    # Deliberately exceed base array sizes: staging and Python state matter too.
    host_per_cell = 8192 + 768 * contacts + 32 * species
    device_per_cell = 1024 + 256 * contacts + 16 * species
    partitioned = getattr(getattr(sim, 'CLDevicePool', None), 'partitioned', False)
    shares = sim.clDeviceWeights if partitioned else [1.0] * len(sim.CLDevices)
    candidates = [int(available * fraction / host_per_cell)]
    for device, share in zip(sim.CLDevices, shares):
        # Fourfold dependency reserve; maximum single buffer also constrains us.
        candidates.append(int(device.global_mem_size * fraction / (4 * device_per_cell * share)))
        candidates.append(int(device.max_mem_alloc_size * fraction / (4 * 32 * contacts * share)))
    candidates.append((2**31 - 1) // (contacts * 8))
    return min(candidates)


def configuration(sim):
    def value(name, default, convert=int):
        return convert(os.environ.get('CM_STRESS_' + name, str(default)))
    cfg = dict(initial_cells=value('INITIAL_CELLS', 1),
               max_contacts=value('MAX_CONTACTS', 32), species=value('SPECIES', 4),
               seed=value('SEED', 12345), growth_rate=value('GROWTH_RATE', 2.0, float),
               report_every=value('REPORT_EVERY', 10),
               memory_fraction=value('MEMORY_FRACTION', 0.5, float))
    if any(cfg[key] < 1 for key in ('initial_cells', 'species', 'report_every')):
        raise ValueError('Initial cells, species and report interval must be positive')
    if cfg['max_contacts'] < 8:
        raise ValueError('max_contacts must be at least 8 (solver scratch requirement)')
    if not 0 < cfg['memory_fraction'] <= 0.8:
        raise ValueError('Memory fraction must be in (0, 0.8]')
    if not math.isfinite(cfg['growth_rate']) or cfg['growth_rate'] <= 0:
        raise ValueError('Growth rate must be positive and finite')
    target = os.environ.get('CM_STRESS_MAX_CELLS', '100000')
    cfg['max_cells'] = (capacity_estimate(sim, cfg['max_contacts'], cfg['species'],
                                        cfg['memory_fraction'], host_available_bytes())
                        if target == 'auto' else int(target))
    if cfg['max_cells'] < cfg['initial_cells']:
        raise ValueError('Capacity is below initial population; reduce initial cells or increase available memory')
    if cfg['max_cells'] * cfg['max_contacts'] * 8 >= 2**31:
        raise ValueError('Capacity exceeds supported signed 32-bit contact/vector indexing')
    cfg['max_sqs'] = value('MAX_SQS', max(192**2, 4 * cfg['max_cells']))
    if not 0 < cfg['max_sqs'] < 2**31:
        raise ValueError('max_sqs must be positive and fit signed 32-bit indexing')
    return cfg


def founder_geometry(count, rng):
    """One founder by default; optional separated, unoriented disk inoculum."""
    if count == 1:
        yield (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)
        return
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(count):
        radius = 6.0 * math.sqrt(index)
        angle = index * golden_angle
        orientation = rng.uniform(0.0, 2.0 * math.pi)
        yield ((radius * math.cos(angle), radius * math.sin(angle), 0.0),
               (math.cos(orientation), math.sin(orientation), 0.0))


def setup(sim):
    from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
    from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator
    from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
    import numpy as np

    global _sim, _cfg, _rng, _previous, _last_report
    _sim, _cfg = sim, configuration(sim)
    _rng = random.Random(_cfg['seed'])
    random.seed(_cfg['seed'])
    np.random.seed(_cfg['seed'])

    class CapacityPhysics(CLBacterium):
        def update_grid(self):
            super().update_grid()
            if self.n_sqs > self.max_sqs:
                raise MemoryError('Spatial grid needs %d bins; raise CM_STRESS_MAX_SQS (currently %d)'
                                  % (self.n_sqs, self.max_sqs))

    print('STRESS configuration:', json.dumps(_cfg, sort_keys=True), flush=True)
    print('Capacity is a target, not measured free VRAM. Population growth stops at the target.', flush=True)
    physics = CapacityPhysics(sim, max_cells=_cfg['max_cells'],
                              max_contacts=_cfg['max_contacts'], max_sqs=_cfg['max_sqs'],
                              jitter_z=False, printing=False)
    integrator = CLEulerIntegrator(sim, _cfg['species'], _cfg['max_cells'])
    sim.init(physics, ModuleRegulator(sim), None, integrator)

    # Match the simple-growth examples: a single cell at the origin grows
    # into a colony. Optional multiple founders are spaced on a disk, with
    # varied orientations, instead of imposing a rectangular aligned lattice.
    for pos, direction in founder_geometry(_cfg['initial_cells'], _rng):
        sim.addCell(cellType=0, pos=pos, dir=direction)
    if sim.is_gui:
        from CellModeller.GUI import Renderers
        sim.addRenderer(Renderers.GLBacteriumRenderer(sim))
    sim.setSaveOutput(False)  # Large pickle output would dominate this test.
    sim.CLWorkStats.clear()
    _previous, _last_report = {}, time.perf_counter()


def init(cell):
    cell.targetVol = 2.5 + _rng.uniform(0.0, 0.5)
    cell.growthRate = _cfg['growth_rate']
    cell.species[:] = 0.2
    cell.color = [0.1, 1.0, 0.3]


def update(cells):
    # Simulator divides all flagged parents in one step. Reserve one additional
    # slot for EACH approved division, not just one margin for the whole step.
    slots = max(0, _cfg['max_cells'] - len(cells))
    for cell in cells.values():
        cell.divideFlag = bool(slots and cell.volume > cell.targetVol)
        if cell.divideFlag:
            slots -= 1
        cell.growthRate = _cfg['growth_rate']
    if not slots:
        for cell in cells.values():
            cell.growthRate = 0.0
    if _sim.is_gui and _sim.stepNum and _sim.stepNum % _cfg['report_every'] == 0:
        report(_sim)


def divide(parent, daughter1, daughter2):
    for daughter in (daughter1, daughter2):
        daughter.targetVol = 2.5 + _rng.uniform(0.0, 0.5)


def specRateCL():
    return '\n'.join('rates[%d] = 0.2f + 0.1f * species[%d] - 0.3f * species[%d];'
                     % (i, (i + 1) % _cfg['species'], i) for i in range(_cfg['species']))


def report(sim, stream=None):
    global _previous, _last_report
    now = time.perf_counter()
    counts = {key: list(values) for key, values in sim.CLWorkStats.items()}
    delta = {key: [int(v - old) for v, old in zip(values, _previous.get(key, [0]*len(values)))]
             for key, values in counts.items()}
    memory = {key: [int(row['bytes']) for row in rows] for key, rows in sim.CLMemoryStats.items()}
    record = dict(step=sim.stepNum, cells=len(sim.cellStates), target=_cfg['max_cells'],
                  interval_seconds=now - _last_report, work_items=delta,
                  last_stage_planned_bytes=memory, host_available_bytes=host_available_bytes())
    print('STRESS step=%d cells=%d/%d interval=%.3fs contacts_work=%s' %
          (record['step'], record['cells'], record['target'], record['interval_seconds'],
           delta.get('physics.find_contacts', 'single-device')), flush=True)
    if stream:
        stream.write(json.dumps(record) + '\n')
        stream.flush()
    _previous, _last_report = counts, now
    return record


def validate_state(sim):
    import numpy as np
    n = len(sim.cellStates)
    if n > sim.phys.max_cells or sim.phys.n_cells != n:
        raise RuntimeError('Population/index capacity invariant failed')
    arrays = (sim.phys.cell_centers[:n].view(np.float32),
              sim.phys.cell_dirs[:n].view(np.float32), sim.phys.cell_lens[:n],
              sim.integ.specLevel[:n])
    if any(not np.isfinite(array).all() for array in arrays):
        raise FloatingPointError('Nonfinite geometry or species detected')
    if np.any(sim.phys.cell_lens[:n] <= 0):
        raise FloatingPointError('Nonpositive cell length detected')



def checkpoint_signature():
    """Reject resume against different simulation source, including dirty edits."""
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    paths = sorted((root / 'CellModeller').rglob('*.py'))
    paths += sorted((root / 'CellModeller').rglob('*.cl'))
    paths.append(Path(__file__).resolve())
    for path in paths:
        digest.update(str(path.relative_to(root)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def save_checkpoint(path, sim, model, reached, hold_steps):
    """Replace one checkpoint atomically, only after a fully completed step."""
    import numpy as np
    data = dict(format='cellmodeller-stress-v1', source=checkpoint_signature(),
                config=dict(model._cfg), dt=sim.dt, reached=reached,
                hold_steps=hold_steps, stepNum=sim.stepNum,
                cellStates=sim.cellStates, lineage=sim.lineage,
                next_id=sim._next_id, next_idx=sim._next_idx,
                specData=sim.integ.levels,
                model_rng=model._rng.getstate(), python_rng=random.getstate(),
                numpy_rng=np.random.get_state())
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.name + '.',
                                         suffix='.tmp', delete=False) as stream:
            tmp = stream.name
            pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
        tmp = None
    finally:
        if tmp is not None:
            os.unlink(tmp)
    print('CHECKPOINT saved: %s (next step %d, %d cells)' %
          (path, sim.stepNum, len(sim.cellStates)), flush=True)


def read_checkpoint(path):
    # Pickle is executable data: only load checkpoints you created/trust.
    with open(path, 'rb') as stream:
        data = pickle.load(stream)
    if data.get('format') != 'cellmodeller-stress-v1':
        raise ValueError('Not a resumable headless stress checkpoint (GUI pickles are unsupported)')
    if data['source'] != checkpoint_signature():
        raise ValueError('Checkpoint source differs: restore the original CellModeller checkout')
    return data


def restore_checkpoint(sim, model, data):
    import numpy as np
    sim.loadFromPickle(data)
    sim._next_id, sim._next_idx = data['next_id'], data['next_idx']
    # The generic pickle loader does not restore this integrator counter.
    sim.integ.nCells = len(sim.cellStates)
    model._rng.setstate(data['model_rng'])
    random.setstate(data['python_rng'])
    np.random.set_state(data['numpy_rng'])
    model._previous, model._last_report = {}, time.perf_counter()
    sim.CLWorkStats.clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--platform', type=int, default=0)
    parser.add_argument('--devices', default='auto', help="'auto' or comma-separated indices; one index tests single GPU")
    parser.add_argument('--weights', help='Comma-separated positive weights; omitted means equal')
    parser.add_argument('--gpu-memory', choices=['partitioned', 'replicated'], default='partitioned')
    parser.add_argument('--max-cells', default='100000', help="Population target, or 'auto' for a memory estimate")
    parser.add_argument('--initial-cells', type=int, default=1, help='Founders; default one cell like the growth examples')
    parser.add_argument('--max-contacts', type=int, default=32)
    parser.add_argument('--max-sqs', type=int, help='Grid capacity; default max(192**2, 4*max_cells)')
    parser.add_argument('--species', type=int, default=4)
    parser.add_argument('--seed', type=int, default=12345)
    parser.add_argument('--growth-rate', type=float, default=2.0)
    parser.add_argument('--memory-fraction', type=float, default=0.5, help='Auto estimate fraction, (0, 0.8]; not a runtime quota')
    parser.add_argument('--report-every', type=int, default=10)
    parser.add_argument('--dt', type=float, default=0.025)
    parser.add_argument('--steps', type=int, default=10000, help='Maximum additional steps in this invocation')
    parser.add_argument('--hold-steps', type=int, default=100, help='Steps at target with growth disabled')
    parser.add_argument('--log', default='multigpu-stress.jsonl', help='Exclusive-create JSONL log path')
    parser.add_argument('--checkpoint', help='Atomic rolling checkpoint filename; enables safe stop')
    parser.add_argument('--checkpoint-every', type=int, default=100, help='Completed steps between checkpoints')
    parser.add_argument('--resume', help='Resume a trusted headless checkpoint; restores workload, dt and hold progress')
    args = parser.parse_args()
    if args.checkpoint_every < 1:
        parser.error('checkpoint-every must be positive')
    resumed = read_checkpoint(args.resume) if args.resume else None
    if resumed:
        for key, value in resumed['config'].items():
            setattr(args, key, value)
        args.dt, args.hold_steps = resumed['dt'], resumed['hold_steps']
        if not args.checkpoint:
            args.checkpoint = args.resume
        print('RESUME: saved workload, dt and hold settings override command-line values', flush=True)
    if args.checkpoint and Path(args.checkpoint).exists() and not (
            args.resume and Path(args.checkpoint).resolve() == Path(args.resume).resolve()):
        parser.error('Checkpoint already exists; use --resume or a new checkpoint filename')
    stop_requested = []
    old_handlers = {}

    if not math.isfinite(args.dt) or args.dt <= 0 or args.steps < 1 or args.hold_steps < 1:
        parser.error('dt, steps and hold-steps must be positive')
    try:
        devices = 'auto' if args.devices == 'auto' else [int(i) for i in args.devices.split(',')]
        weights = None if args.weights is None else [float(i) for i in args.weights.split(',')]
    except ValueError:
        parser.error('Invalid device indices or weights')
    for key in ('max_cells', 'initial_cells', 'max_contacts', 'max_sqs', 'species', 'seed',
                'growth_rate', 'memory_fraction', 'report_every'):
        if getattr(args, key) is not None:
            os.environ['CM_STRESS_' + key.upper()] = str(getattr(args, key))
    # Prefer this checkout when executing the example directly.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from CellModeller.Simulator import Simulator
    with open(args.log, 'x') as stream:
        sim = None
        try:
            if args.checkpoint:
                def request_stop(signum, frame):
                    stop_requested.append(signum)
                for signum in (signal.SIGINT, signal.SIGTERM):
                    old_handlers[signum] = signal.signal(signum, request_stop)
            sim = Simulator(str(Path(__file__).resolve()), args.dt, saveOutput=False,
                            clPlatformNum=args.platform, clDeviceNums=devices,
                            clDeviceWeights=weights, clMultiGPUMemory=args.gpu_memory)
            model = sim.module  # Simulator imports a separate model module.
            stream.write(json.dumps(dict(event='start', config=model._cfg, arguments=vars(args),
                                         devices=[dict(name=d.name, memory=int(d.global_mem_size),
                                                       driver=d.driver_version) for d in sim.CLDevices],
                                         weights=sim.clDeviceWeights)) + '\n')
            stream.flush()
            reached = resumed['reached'] if resumed else None
            if resumed:
                restore_checkpoint(sim, model, resumed)
                model.validate_state(sim)
                stream.write(json.dumps(dict(event='resume', step=sim.stepNum,
                                             checkpoint=args.resume)) + '\n')
                stream.flush()
            if args.checkpoint:
                save_checkpoint(args.checkpoint, sim, model, reached, args.hold_steps)
            stopped = False
            for _ in range(args.steps):
                if stop_requested:
                    stopped = True
                    break
                if reached is not None and sim.stepNum - reached >= args.hold_steps:
                    break
                sim.step()
                if sim.stepNum % model._cfg['report_every'] == 0:
                    model.validate_state(sim)
                    model.report(sim, stream)
                if len(sim.cellStates) == model._cfg['max_cells']:
                    if reached is None:
                        reached = sim.stepNum
                if args.checkpoint and sim.stepNum % args.checkpoint_every == 0:
                    save_checkpoint(args.checkpoint, sim, model, reached, args.hold_steps)
                if stop_requested:
                    stopped = True
                    break
            model.validate_state(sim)
            model.report(sim, stream)
            if args.checkpoint:
                save_checkpoint(args.checkpoint, sim, model, reached, args.hold_steps)
            status = 'stopped_checkpointed' if stopped else 'target_held' if reached is not None and sim.stepNum - reached >= args.hold_steps else 'step_limit'
            stream.write(json.dumps(dict(event='end', status=status, cells=len(sim.cellStates))) + '\n')
            print('STRESS finished:', status, '(step_limit does not establish maximum capacity)', flush=True)
        except (Exception, KeyboardInterrupt) as error:
            stream.write(json.dumps(dict(event='failure', error=repr(error),
                                         step=getattr(sim, 'stepNum', None),
                                         cells=len(sim.cellStates) if sim else None)) + '\n')
            stream.flush()
            raise
        finally:
            for signum, handler in old_handlers.items():
                signal.signal(signum, handler)


if __name__ == '__main__':
    main()
