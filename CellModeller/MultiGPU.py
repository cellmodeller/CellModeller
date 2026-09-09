"""Explicit, synchronized work partitioning for simulation OpenCL kernels.

Each worker has private buffers. Only its owned output rows are merged, so
contacts across partition boundaries still use the complete neighbor geometry.
"""
from collections import OrderedDict
import math
import numbers

import pyopencl as cl


def device_indices(devices, selection, default=0):
    if selection is None:
        indices = [default]
    elif isinstance(selection, str):
        if selection != 'auto':
            raise ValueError("clDeviceNums must be a sequence of indices or 'auto'")
        indices = [i for i, d in enumerate(devices) if d.type & cl.device_type.GPU]
        if not indices:
            indices = [default]
    else:
        try:
            indices = list(selection)
        except TypeError as error:
            raise ValueError("clDeviceNums must be a sequence of indices or 'auto'") from error
    if not indices or any(isinstance(i, bool) or not isinstance(i, numbers.Integral)
                          or i < 0 or i >= len(devices) for i in indices):
        raise ValueError('OpenCL device indices must be a nonempty list of valid indices')
    if len(set(indices)) != len(indices):
        raise ValueError('OpenCL device indices must be unique')
    if len(indices) > 1 and any(not devices[i].type & cl.device_type.GPU for i in indices):
        raise ValueError('Multi-GPU execution requires GPU devices on the selected platform')
    return indices


def normalized_weights(weights, count):
    weights = [1.0] * count if weights is None else list(weights)
    if len(weights) != count:
        raise ValueError('clDeviceWeights must have one entry per selected device')
    weights = [float(w) for w in weights]
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError('clDeviceWeights must be finite and positive')
    scale = max(weights)
    weights = [w / scale for w in weights]
    total = sum(weights)
    return [w / total for w in weights]


def partitions(count, weights):
    """Contiguous ranges, including empty ranges, covering every cell once."""
    # Largest-remainder apportionment avoids cumulative floating-point error
    # (e.g. six equal GPUs receiving 99/101 rows instead of 100/100).
    quotas = [count*weight for weight in weights]
    sizes = [int(quota) for quota in quotas]
    remaining = count-sum(sizes)
    order = sorted(range(len(weights)), key=lambda i: (-(quotas[i]-sizes[i]), i))
    for index in order[:remaining]:
        sizes[index] += 1
    start = 0
    for size in sizes:
        yield start, size
        start += size


class DevicePool:
    """Bounded replica cache shared by all stages of one simulation.

    Primary buffers remain authoritative between calls (including user .set()
    calls). Cached replicas are refreshed before every launch. This is allocation
    reuse, not an assumption that an external writer left the buffer unchanged.
    """
    def __init__(self, queues, weights, min_cells, cache_bytes=256*1024*1024):
        self.queues = queues
        self.weights = weights
        self.min_cells = min_cells
        self.cache_bytes = cache_bytes
        self.caches = [OrderedDict() for _ in queues]
        self.cache_sizes = [0 for _ in queues]
        self.stats = {}

    def replica(self, device, original):
        key = getattr(original, 'int_ptr', id(original))
        cache = self.caches[device]
        if key in cache:
            cache.move_to_end(key)
            return cache[key][1]
        replica = cl.Buffer(self.queues[device].context, cl.mem_flags.READ_WRITE, original.size)
        if original.size <= self.cache_bytes:
            while cache and self.cache_sizes[device]+original.size > self.cache_bytes:
                _, (old, _) = cache.popitem(last=False)
                self.cache_sizes[device] -= old.size
            # Retaining the original prevents reuse of its pointer as another
            # allocation while this cache entry exists. The cache is bounded.
            cache[key] = (original, replica)
            self.cache_sizes[device] += original.size
        return replica

    def record(self, name, ranges, trailing=1):
        counts = self.stats.setdefault(name, [0]*len(self.queues))
        for i, (_, size) in enumerate(ranges):
            counts[i] += size*trailing

    def clear_cache(self):
        for queue in self.queues:
            queue.finish()
        for cache in self.caches:
            cache.clear()
        self.cache_sizes[:] = [0]*len(self.queues)


class DistributedProgram:
    """Split dimension zero; each output specifies (argument, row bytes, offset).

    Layouts must be audited: each work item may write only its owned output row.
    Kernels may read arbitrary neighbors from the full input snapshot. Reduction
    kernels are handled separately by DistributedReduction below.
    """
    def __init__(self, program, pool, layouts, label=''):
        self.program = program
        self.pool = pool
        self.output_layouts = layouts
        self.label = label
        self.kernels = {}

    def __getattr__(self, name):
        if name in self.kernels:
            return self.kernels[name]
        kernel = getattr(self.program, name)
        if name not in self.output_layouts:
            return kernel

        def run(queue, global_size, local_size, *args, **kwargs):
            pool = self.pool
            count = global_size[0]
            if any(size == 0 for size in global_size):
                return cl.enqueue_marker(queue)
            if len(pool.queues) == 1 or count < pool.min_cells:
                pool.record(self.label+name, [(0, count)]+[(0, 0)]*(len(pool.queues)-1), math.prod(global_size[1:]))
                return kernel(queue, global_size, local_size, *args, **kwargs)
            if queue != pool.queues[0] or local_size is not None or kwargs:
                raise ValueError('Distributed kernels require the primary queue and an unoffset launch')
            outputs = self.output_layouts[name](args)
            for index, width, offset in outputs:
                if width < 0 or offset < 0 or offset+count*width > args[index].size:
                    raise ValueError('Distributed output layout exceeds buffer size')
            ranges = list(partitions(count, pool.weights))
            # Complete all input snapshots before ANY kernel writes, including
            # the primary. Private replicas also preserve aliased arguments.
            queue.finish()
            workers = []
            for device, (start, size) in enumerate(ranges[1:], 1):
                if not size:
                    continue
                worker_args = list(args)
                copied = {}
                for index, arg in enumerate(args):
                    if isinstance(arg, cl.Buffer):
                        key = getattr(arg, 'int_ptr', id(arg))
                        if key not in copied:
                            copied[key] = pool.replica(device, arg)
                            cl.enqueue_copy(queue, copied[key], arg)
                        worker_args[index] = copied[key]
                workers.append((pool.queues[device], start, size, worker_args))
            queue.finish()
            events = []
            start, size = ranges[0]
            if size:
                events.append(kernel(queue, (size,)+tuple(global_size[1:]), None, *args,
                                     global_offset=(start,)+(0,)*(len(global_size)-1)))
                queue.flush()
            for worker_queue, start, size, worker_args in workers:
                events.append(kernel(worker_queue, (size,)+tuple(global_size[1:]), None,
                                     *worker_args, global_offset=(start,)+(0,)*(len(global_size)-1)))
                worker_queue.flush()
            for event in events:
                event.wait()
            for worker_queue, start, size, worker_args in workers:
                for index, width, offset in outputs:
                    if width:
                        cl.enqueue_copy(queue, args[index], worker_args[index],
                                        src_offset=offset+start*width, dst_offset=offset+start*width,
                                        byte_count=size*width)
            event = cl.enqueue_marker(queue)
            event.wait()  # Replicas stay alive until all merges finish.
            pool.record(self.label+name, ranges, math.prod(global_size[1:]))
            return event
        self.kernels[name] = run
        return run


def rows(*specs):
    """Output (argument, bytes) pairs, with a callable for dynamic strides."""
    return lambda args: [(i, width(args) if callable(width) else width, 0) for i, width in specs]


def contact_outputs(first, contact_arg, overlap):
    def layout(args):
        stride = int(args[contact_arg])
        widths = [4] + [stride*b for b in (4, 4, 4, 16, 16, 4, 4)]
        if overlap:
            widths.append(stride*4)
        return [(i, width, 0) for i, width in enumerate(widths, first)]
    return layout


CONTACT_LAYOUTS = {'find_plane_contacts': (10, 1, False),
                   'find_sphere_contacts': (11, 1, False),
                   'find_contacts': (15, 7, True)}
PHYSICS_LAYOUTS = {name: contact_outputs(*layout) for name, layout in CONTACT_LAYOUTS.items()}
PHYSICS_LAYOUTS.update({
    'bin_cells': rows((6, 4)),
    'collect_tos': rows((14, lambda a: int(a[7])*4), (15, 4)),
    'build_matrix': rows((10, lambda a: int(a[0])*32), (11, lambda a: int(a[0])*32)),
    'calculate_Bx': rows((6, lambda a: int(a[0])*4)),
    'calculate_BTBx': rows((7, 32)),
    'calculate_Mx': rows((6, 32)),
    'calculate_Minv_x': rows((6, 32)),
    'predict': rows((6, 16), (7, 16), (8, 4)),
    'integrate': rows((0, 16), (1, 16), (2, 4), (3, 16), (4, 16), (5, 4)),
    'add_impulse': rows((6, 16), (7, 16), (9, 4)),
})
EULER_LAYOUTS = {
    'speciesRates': rows((7, lambda a: int(a[0])*4)),
    'diluteSpecs': rows((3, lambda a: int(a[0])*4)),
}
SIGNAL_LAYOUTS = {
    'gridCells': rows((10, 32), (11, 32)),
    'setCellSignals': rows((8, lambda a: int(a[0])*4)),
    # The user species function has a writable species pointer in these
    # integrators, so preserve changes to species as well as to rates.
    'speciesRates': rows((6, lambda a: int(a[1])*4), (8, lambda a: int(a[1])*4)),
    'signalRates': rows((9, lambda a: int(a[0])*32)),
    'diluteSpecs': rows((3, lambda a: int(a[0])*4)),
    'speciesDT': rows((1, lambda a: int(a[0])*4)),
    'setCellSigImplicit': rows((4, lambda a: int(a[5])*8)),
    'setCellSignalsImplicit': rows((9, lambda a: int(a[0])*8)),
    'speciesRatesImplicit': rows((6, lambda a: int(a[1])*4), (8, lambda a: int(a[1])*4)),
    'signalRatesImplicit': rows((9, lambda a: int(a[0])*32)),
}


class ContactProgram(DistributedProgram):
    """Compatibility adapter for callers opting into contact-only dispatch."""
    layouts = CONTACT_LAYOUTS

    def __init__(self, program, queues, weights, min_cells):
        super().__init__(program, DevicePool(queues, weights, min_cells),
                         {name: contact_outputs(*layout) for name, layout in self.layouts.items()})


def distribute_program(sim, program, layouts, label, source=None):
    pool = getattr(sim, 'CLDevicePool', None)
    if getattr(pool, 'partitioned', False):
        from .PartitionedGPU import PartitionedProgram
        return PartitionedProgram(source if source is not None else program.source, pool, layouts, label)
    return DistributedProgram(program, pool, layouts, label) if pool is not None else program


class DistributedElementwise:
    """Adapter for the row-local ElementwiseKernels used by bacterial physics.

    Explicit per-array offsets support views, including nonzero starts, without
    relying on device-specific sub-buffer alignment. Output is the first array.
    """
    def __init__(self, original, pool, arguments, operation, name):
        import re
        self.original = original
        self.pool = pool
        self.name = name
        self.arguments = []
        declarations, aliases = [], []
        for declaration in arguments.split(','):
            match = re.fullmatch(r'\s*(?:__global\s+)?(?:const\s+)?(\w+)\s*(\*?)\s*(\w+)\s*', declaration)
            if not match:
                raise ValueError('Unsupported elementwise argument: '+declaration)
            dtype, pointer, argname = match.groups()
            self.arguments.append((dtype, bool(pointer), argname))
            if pointer:
                declarations.extend(['__global %s *%s_data' % (dtype, argname),
                                     'ulong %s_offset' % argname])
                aliases.append('__global %s *%s = %s_data + %s_offset;' %
                               (dtype, argname, argname, argname))
            else:
                declarations.append('%s %s' % (dtype, argname))
        source = '__kernel void %s(%s) { size_t i=get_global_id(0); %s %s; }' % (
            name, ', '.join(declarations), ' '.join(aliases), operation)
        self.program = cl.Program(pool.queues[0].context, source).build(cache_dir=False)

    def __call__(self, *args, **kwargs):
        import numpy as np
        output = args[0]
        count = output.size
        if len(self.pool.queues) == 1 or count < self.pool.min_cells:
            self.pool.record(self.name, [(0, count)]+[(0, 0)]*(len(self.pool.queues)-1))
            return self.original(*args, **kwargs)
        if kwargs or len(args) != len(self.arguments):
            raise ValueError('Distributed elementwise launch expects positional array/scalar arguments')
        raw = []
        for arg, (dtype, pointer, _) in zip(args, self.arguments):
            if pointer:
                itemsize = {'float': 4, 'int': 4, 'float4': 16, 'float8': 32}[dtype]
                if (not arg.flags.c_contiguous or arg.nbytes < count*itemsize
                        or arg.offset % itemsize):
                    raise ValueError('Distributed elementwise arrays must be contiguous and cover the output')
                raw.extend([arg.base_data, np.uint64(arg.offset//itemsize)])
            else:
                raw.append({'float': np.float32, 'int': np.int32}[dtype](arg))
        width = output.dtype.itemsize
        # Bind per-call output offset; the generated code does the same for
        # inputs. A fresh lightweight wrapper reuses the compiled program/pool.
        layout = lambda values: [(0, width, int(values[1])*width)]
        program = DistributedProgram(self.program, self.pool, {self.name: layout})
        return getattr(program, self.name)(self.pool.queues[0], (count,), None, *raw)


class DistributedReduction:
    """Run the existing associative reduction on every partition; sum scalars.

    The CG convergence decision is still made once on the host, from the global
    sum. No GPU computes a duplicate full-vector reduction.
    """
    def __init__(self, original, pool, dtype, name):
        self.original = original
        self.pool = pool
        self.dtype = dtype
        self.name = name

    def __call__(self, *args, **kwargs):
        import numpy as np
        import pyopencl.array as cl_array
        count = args[0].size
        pool = self.pool
        if len(pool.queues) == 1 or count < pool.min_cells:
            pool.record(self.name, [(0, count)]+[(0, 0)]*(len(pool.queues)-1))
            return self.original(*args, **kwargs)
        if kwargs or any(not a.flags.c_contiguous or a.size != count for a in args):
            raise ValueError('Distributed reductions require equal-sized contiguous arrays')
        ranges = list(partitions(count, pool.weights))
        pool.queues[0].finish()
        workers = []
        for device, (start, size) in enumerate(ranges):
            if not size:
                continue
            views = []
            copied = {}
            for arg in args:
                buffer = arg.base_data
                if device:
                    key = getattr(buffer, 'int_ptr', id(buffer))
                    if key not in copied:
                        copied[key] = pool.replica(device, buffer)
                        cl.enqueue_copy(pool.queues[0], copied[key], buffer)
                    buffer = copied[key]
                views.append(cl_array.Array(pool.queues[device], (size,), arg.dtype,
                                            data=buffer, offset=arg.offset+start*arg.dtype.itemsize))
            workers.append((pool.queues[device], views))
        pool.queues[0].finish()
        partials = []
        for queue, views in workers:
            partials.append(self.original(*views, queue=queue))
            queue.flush()
        # Accumulate float partials in double precision before casting back to
        # the legacy dtype; integer sums use Python integers to avoid overflow.
        values = [partial.get().item() for partial in partials]
        total = math.fsum(values) if np.dtype(self.dtype).kind == 'f' else sum(values)
        result = cl_array.to_device(pool.queues[0], np.asarray(total, dtype=self.dtype))
        pool.record(self.name, ranges)
        return result
