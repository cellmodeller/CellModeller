# Multi-GPU simulations and larger populations

Multi-GPU simulations now default to **partitioned GPU memory**. The host keeps
canonical cell, contact, solver, and integration arrays. Each GPU processes its
owned cells and receives only the neighboring data required by that stage.
There is no full-state allocation on GPU 0 and no full-state replica on every
GPU in this mode.

This removes the previous requirement that the complete GPU working set fit
on one device. Larger populations can fit across multiple GPUs, subject to
host RAM, each partition's working set, and the size of its neighbor data.
This is host-staged distributed computation, not a transparent unified VRAM
address space or a promise of linear capacity/speed scaling.

## Enable partitioned memory

```python
from CellModeller.Simulator import Simulator

sim = Simulator(
    'my_large_model.py', 0.025,
    clPlatformNum=0,
    clDeviceNums=[0, 1],       # Or 'auto' for all GPUs on this platform
    clDeviceWeights=None,     # Equal shares, appropriate for identical GPUs
    clMultiGPUMemory='partitioned',  # Default for multi-GPU simulations
)
sim.step()
print(sim.CLWorkStats)
print(sim.CLMemoryStats)
```

Increase the model's total capacity too. For example, in its `setup(sim)`:

```python
biophys = CLBacterium(sim, max_cells=200000, max_contacts=24, max_sqs=512**2)
```

Set each integrator's `maxCells` to at least the same total capacity. These
values remain model limits; selecting more GPUs does not silently change them.
Large colonies can also require a larger `max_sqs` spatial grid. Contact density
must fit `max_contacts`; partitioned contact kernels detect capacity overflow
and abort the stage rather than knowingly commit invalid contact data.

The built-in batch runner now accepts a larger stopping limit:

```bash
python Scripts/batch.py my_large_model.py 0 --devices 0,1 --gpu-memory partitioned --max-cells 200000
```

The runner stops 256 cells below that limit, preserving its existing division
buffer, and checks that the model's `max_cells` is large enough.

## What is stored on each GPU

Spatial bisection divides active cells according to the selected workload
weights. Equal weights assign equal counts, within one cell. Ownership is
recomputed at each physics substep as cells move or divide. Global cell IDs
remain unchanged, so migration does not alter lineage or regulator state.

- Cell-local operations upload only owned rows.
- Contact detection uploads owned contact tables and the cell geometry in the
  neighboring grid bins visited by those cells. Spatial indices are sparse too.
- Matrix construction/products gather only the cell-vector entries referenced
  by the partition's contacts.
- Transpose products gather the precise incoming contact entries they need;
  they do not replicate the contact tables of every neighboring cell.
- Dot products and contact-count reductions upload local input slices and
  combine small partial sums on the host.
- Species buffers are sliced by owned cell. Signal interpolation uploads only
  the referenced grid nodes for each signal, not the full diffusion grid.

Each device has its own OpenCL context in partitioned mode. This prevents
implicit driver placement of all workers' allocations on a single primary
device. All exchanges are staged through host RAM. Kernels remap global indices
into compact local buffers. Missing neighbor entries or contact-capacity errors
abort the stage. Owned results are committed only after every worker succeeds.
GPU scratch buffers are released when the stage completes; full canonical
state resides in host RAM between stages.

## Covered computation

Partitioned execution covers all built-in `CLBacterium` GPU program kernels,
its conjugate-gradient vector operations and reductions, geometry calculations,
`CLFixedPosition` volume growth, and the OpenCL kernels in `CLEulerIntegrator`,
`CLEulerSigIntegrator`, and `CLCrankNicIntegrator`.

Python regulation, cell-state bookkeeping, NumPy sorting/aggregation, and SciPy
diffusion/convolution remain host operations. Their memory requirements still
scale with the total population/grid. Host canonical arrays and temporary
staging copies require additional RAM. GPU transfer and index-remapping overhead
can make this mode slower even though it supports a larger working set.

Custom code that passes built-in `*_dev.data` directly to an arbitrary OpenCL
kernel must be adapted to the partitioned interface. These handles refer to
host-backed arrays in partitioned mode; they do not silently allocate a full
GPU buffer. Rate functions must write only their current cell's species/rates.
Keep `clMultiGPUMemory='replicated'` for older custom GPU code while adapting it.

## Preferences and compatibility

- `clDeviceNums=None`: existing single-device behavior, using `clDeviceNum`.
- `clDeviceNums=[0, 1]`: select unique GPU indices on `clPlatformNum`.
- `clDeviceNums='auto'`: all GPUs on that platform. With one GPU, the normal
  single-GPU path is used. With no GPUs, the existing `clDeviceNum` fallback is used.
- `clDeviceWeights=None`: equal counts. Positive finite weights such as `[2, 1]`
  assign approximately two thirds of the rows to the first GPU. Weights follow
  selection order. Equal counts do not guarantee identical utilization because
  contact complexity, device clocks, and transfers vary.
- `clMultiGPUMemory='partitioned'`: default when more than one GPU is selected.
- `clMultiGPUMemory='replicated'`: previous compute-distribution backend; full
  buffers still need to fit on a device. This mode retains a bounded replica cache.
- `clMultiGPUMinCells=1024`: applies to **replicated mode only**. Partitioned
  mode keeps work split over the selected GPUs, even for small launches, to
  avoid reintroducing a GPU 0 capacity bottleneck. Empty partitions are skipped.

Both modes currently select GPUs on a single OpenCL platform. A single CPU
OpenCL device remains supported by the normal single-device path.

The GUI exposes memory mode, GPU selection, weights, and threshold. Settings
survive viewer reset/pickle loading, but are not stored in simulation pickles
or persisted across application restarts. Batch equivalents are `--devices`,
`--device-weights`, `--gpu-memory`, and `--multi-gpu-min-cells`.

## Inspect distribution and memory

`sim.CLWorkStats` records cumulative work items per stage in device-selection
order. A 2D matrix kernel includes its contact dimension. Clear this dictionary
after initialization to measure a particular interval. Counts are not timings
or FLOP measurements.

`sim.CLMemoryStats` in partitioned mode records the most recent launch per
stage: owned row counts and planned scratch bytes per GPU. Program kernels also
report each buffer's row count, row width, and data bytes. Reductions include a
conservative allowance for intermediate partial sums. These values exclude
compiled-program/driver overhead and are not measured whole-device VRAM usage.

Before launching, the backend checks each allocation against the device's
maximum allocation size and the estimated working set against 80% of its
reported global memory. Oversized partitions raise `MemoryError` with the GPU
and requested bytes. Existing external GPU memory usage can still cause an
allocation failure below this bound. These limits use the
[OpenCL device properties](https://documen.tician.de/pyopencl/runtime_platform.html).

Dense neighborhoods or long-range contact graphs can make a partition's halo
large. Full neighbor replication in such a pathological geometry can still be
necessary for an individual input, even though owned contact/solver/output
storage remains partitioned. Two 8 GB GPUs therefore do not guarantee exactly
twice the supported population of one 8 GB GPU. No maximum population or speedup
has been measured on hardware in the implementation environment.

## Validation

### Growing-colony capacity example

For a dedicated graphical setup with no environment variables, launch:

```bash
python Examples/multigpu_stress_gui.py
```

Choose an explicit population target or memory-based estimate, founders,
contacts, species, growth rate, seed and reporting interval. Then select GPUs,
weights and partitioned memory in the normal device dialog. Click **Run** in
the colony window. A separate status window shows population, selected devices,
weights, cumulative contact work and latest planned contact scratch memory.
Status refreshes between simulation steps, so expensive steps may delay it.
The launcher is standalone; use `multigpu_stress.py` for the **Load Model** action.
To change the launcher's settings, close it and launch it again.

Load `Examples/multigpu_stress.py` in the GUI and select your GPUs with
partitioned memory. The default grows 256 founders toward 100,000 cells with
mechanics, contact/CG solver work, division and four intracellular species.
It reserves division slots so a single step cannot exceed the array capacity.
Growth stops at the target; the GUI can continue running to inspect the colony.
Signal diffusion is not part of this example.

For capacity testing, run without rendering or pickle output:

```bash
python Examples/multigpu_stress.py --devices auto --max-cells auto --log stress-auto.jsonl
```

`auto` estimates a target using available host RAM (including exposed Linux
container limits), device memory, allocation limits, weights and a halo reserve.
It is a planning heuristic, not a measurement of free GPU memory or the maximum
supported population. It does not override the backend's allocation guards.
The default estimate uses a 0.5 memory fraction; `--memory-fraction` changes that
estimate up to 0.8. Explicit targets bypass the estimate and let you approach the
practical limit across separate runs. Large founder counts are expensive to
initialize, so increase `--initial-cells` cautiously.

For controlled comparisons use the same explicit population target, initial
population, seed, timestep and species count, with a different log per run:

```bash
python Examples/multigpu_stress.py --devices 0 --max-cells 200000 --log stress-one.jsonl
python Examples/multigpu_stress.py --devices 0,1 --max-cells 200000 --log stress-two.jsonl
```

The runner logs interval times, population, per-stage work-item counts and
planned per-device scratch. It checks finite geometry/species and positive
lengths at reporting intervals and completion. It holds the target population
for 100 steps by default, with a 10,000-step overall bound. `target_held` means
that target completed; `step_limit` does not mean capacity was reached. Failures
are logged and re-raised, rather than silently reducing the workload. Logs are
created exclusively so earlier results are not overwritten. Use external GPU
telemetry alongside these logs; work counts are not utilization measurements.
Full compute utilization or full VRAM occupancy cannot be guaranteed by a model
file, especially with the host-staged partitioned backend.

For an estimated-capacity GUI run, set the environment before launch:

```bash
CM_STRESS_MAX_CELLS=auto python Scripts/CellModellerGUI.py
```

Then load the example and select the GPUs. `--help` documents all CLI options
and their GUI environment equivalents. GUI diagnostics describe the preceding
completed steps. Use the hardware comparison suite below to test numerical
agreement; successful stress execution alone does not establish equivalence.

```bash
python -m unittest discover -s tests -v
```

The scheduling suite runs without GPU dependencies. Partitioned-memory tests
require NumPy and use a fake device runtime. They check that model allocation
uses no GPU buffers, per-device scratch requirements fall as GPUs are added,
spatial ownership preserves IDs, array slices remain correct, halo sets are
sparse, reductions upload local slices, and errors cannot partially commit a
stage. Kernel source rewriting and capacity guards are also checked.

With the normal CellModeller dependencies and two GPUs, the hardware suite
compares real physics and integration kernels with both multi-GPU backends,
and compares complete single/partitioned simulations through solver work,
species integration, and cell division. Hardware tests are included but have
not run here. Validate numerical tolerances, memory use, and performance on
representative colony sizes before relying on a larger production run.
