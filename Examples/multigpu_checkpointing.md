# Resumable headless stress runs

`multigpu_stress.py` supports `--checkpoint FILE`, `--checkpoint-every N`
(default 100 completed steps), and `--resume FILE`. This applies to the headless
stress example, not the GUI example's Save Pickles button.

A checkpoint contains cell states, lineage, species levels, next-cell identifiers,
completed-step position, resolved population/grid configuration, timestep, target
hold progress, and model/Python/NumPy random states. GPU scratch/solver buffers
are rebuilt when loading. This supports continuation, not a guarantee of
bit-for-bit numerical equivalence to an uninterrupted GPU run.

The resolved capacity is restored, including when the original target was `auto`.
Saved workload, timestep and hold settings override their command-line values on
resume. Devices, weights and memory mode are chosen again; use the same hardware
and settings first. The chosen hardware must fit the saved capacity. Resume rejects
changed Python/OpenCL simulation source. Keep the exact checkout and dependencies.
Only load trusted checkpoints: Python pickle can execute code.

## Short hardware smoke test first

From `/opt/CellModeller`, after installing this updated branch:

```bash
SMOKE_DIR=$(mktemp -d "$HOME/cm-resume-smoke-XXXXXXXX")
python -u Examples/multigpu_stress.py --devices 0,1,2,3 \
  --max-cells 32 --max-sqs 36864 --steps 3 \
  --checkpoint "$SMOKE_DIR/latest.pickle" --checkpoint-every 1 \
  --log "$SMOKE_DIR/first.jsonl"
python -u Examples/multigpu_stress.py --devices 0,1,2,3 \
  --resume "$SMOKE_DIR/latest.pickle" --steps 3 \
  --log "$SMOKE_DIR/resumed.jsonl"
```

Expect the first checkpoint to finish at next step 3 and the resumed one at next
step 6, without errors. `step_limit` is expected for this smoke test. It checks
loading and continued stepping, not scientific equivalence.

## Four-device 100,000-cell run

Use the same terminal to create an output directory, then start telemetry:

```bash
cd /opt/CellModeller
RUN_DIR=$(mktemp -d "$HOME/cm-four-gpu-resumable-XXXXXXXX")
echo "$RUN_DIR"
git rev-parse HEAD > "$RUN_DIR/commit.txt"
git diff HEAD > "$RUN_DIR/source-changes.patch"
python -m pip freeze > "$RUN_DIR/python-packages.txt"
nvidia-smi -q > "$RUN_DIR/nvidia-details.txt"
nvidia-smi --query-gpu=timestamp,index,utilization.gpu,utilization.memory,memory.used,power.draw \
  --format=csv -l 1 > "$RUN_DIR/telemetry-01.csv" &
TELEMETRY_PID=$!
python -u Examples/multigpu_stress.py \
  --platform 0 --devices 0,1,2,3 --weights 1,1,1,1 \
  --gpu-memory partitioned --max-cells 100000 --initial-cells 1 \
  --max-contacts 32 --max-sqs 400000 --species 4 --seed 12345 \
  --growth-rate 2 --memory-fraction 0.5 --report-every 10 \
  --dt 0.025 --steps 10000 --hold-steps 100 \
  --checkpoint "$RUN_DIR/latest.pickle" --checkpoint-every 100 \
  --log "$RUN_DIR/stress-01.jsonl" > "$RUN_DIR/console-01.log" 2>&1
# After Python returns (including a completed graceful stop):
kill "$TELEMETRY_PID" 2>/dev/null || true
wait "$TELEMETRY_PID" 2>/dev/null || true
```

In another terminal, use `tail -f /full/path/to/console-01.log` to watch progress.
Keep the launch terminal open. Checkpoint writes add overhead to reported timing;
these timings should not be presented as an identical no-checkpoint benchmark.

## Stop and resume

Press Ctrl+C once in the terminal running Python. With checkpointing enabled,
SIGINT and SIGTERM request a graceful stop. It can take a long step to finish;
wait for `CHECKPOINT saved` and `STRESS finished: stopped_checkpointed` in the
console log and for Python to exit. Signals do not interrupt the solver mid-step.
Do not use `kill -9` or end the hosted session to request a fresh checkpoint.
On a crash/forced termination only the last successfully written checkpoint is
available. Ordinary exceptions do not save potentially half-updated state.

The checkpoint is written to a temporary file, flushed and fsynced, then atomically
replaces the previous checkpoint on the same filesystem. A failed write leaves the
previous file intact. It is a single rolling checkpoint, not a history; back it up
after stopping. Store the run folder in persistent storage before ending a session.
Filesystem persistence across an AxonOS session is separate from this save logic.

On return, set `RUN_DIR` to the existing directory and start a **new** telemetry
file using the same `nvidia-smi` command, then:

```bash
cd /opt/CellModeller
RUN_DIR=/full/path/to/your/existing/run-folder
python -u Examples/multigpu_stress.py \
  --platform 0 --devices 0,1,2,3 --weights 1,1,1,1 \
  --gpu-memory partitioned --resume "$RUN_DIR/latest.pickle" \
  --checkpoint-every 100 --steps 10000 \
  --log "$RUN_DIR/stress-02.jsonl" > "$RUN_DIR/console-02.log" 2>&1
```

Resume updates the same rolling checkpoint unless a new `--checkpoint` path is
given. Logs must have new filenames. `--steps` counts additional steps in this
invocation, and hold progress carries over. Each invocation has its own timing and
work counters; do not combine step deltas as though a resumed log starts at zero.
The start record includes `arguments.resume`; the resume event records its step.
GUI pickles are deliberately rejected because they lack this recovery metadata.
