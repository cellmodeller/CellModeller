"""Host-staged spatial partitions: no full simulation buffer lives on any GPU.

Host state is authoritative between stages. Each launch uploads only owned
rows and its read halo, then commits owned outputs after every device succeeds.
This trades PCIe traffic/host RAM for aggregate GPU capacity. It is deliberately
separate from the backwards-compatible replicated backend.
"""
import math
import re
from dataclasses import dataclass

import numpy as np
import pyopencl as cl

from .MultiGPU import partitions


class Completed:
    def wait(self):
        return None


class HostBuffer:
    def __init__(self, array):
        if not array.flags.c_contiguous:
            raise ValueError('Partitioned kernel buffers must be contiguous')
        self.array = array
        self.size = array.nbytes

    def bytes(self):
        return self.array.reshape(-1).view(np.uint8).reshape(-1)


class PartitionedArray:
    """Subset of the device-array API used by the built-in simulation models.

    Allocation reserves host RAM only. Exposing .data returns a typed host
    handle, never an implicit full-sized GPU allocation. Unsupported custom
    OpenCL calls therefore fail explicitly instead of exhausting GPU 0.
    """
    def __init__(self, array, queue=None):
        self.array = array
        self.queue = queue

    @property
    def data(self):
        return HostBuffer(self.array)

    @property
    def base_data(self):
        return self.data

    @property
    def offset(self):
        return 0

    def __getattr__(self, name):
        if name in ('shape', 'size', 'dtype', 'nbytes', 'flags', 'ndim', 'strides'):
            return getattr(self.array, name)
        raise AttributeError(name)

    def __getitem__(self, key):
        value = self.array[key]
        # Slice views are the normal path; scalar reads retain .get() behavior.
        return PartitionedArray(np.asarray(value), self.queue)

    def __setitem__(self, key, value):
        self.array[key] = value.array if isinstance(value, PartitionedArray) else value

    def get(self, **kwargs):
        if kwargs:
            raise ValueError('Partitioned .get() takes no queue/async options')
        return self.array.copy()

    def set(self, value, **kwargs):
        if kwargs:
            raise ValueError('Partitioned .set() takes no queue/async options')
        np.copyto(self.array, np.asarray(value), casting='same_kind')

    def __mul__(self, value):
        return PartitionedArray(self.array*value, self.queue)

    def __iadd__(self, value):
        self.array += value.array if isinstance(value, PartitionedArray) else value
        return self


class PartitionedArrays:
    @staticmethod
    def zeros(queue, shape, dtype):
        return PartitionedArray(np.zeros(shape, dtype), queue)


def array_module(sim):
    if getattr(getattr(sim, 'CLDevicePool', None), 'partitioned', False):
        return PartitionedArrays
    import pyopencl.array as arrays
    return arrays


def copy_array(queue, destination, source):
    if isinstance(destination, HostBuffer) and isinstance(source, HostBuffer):
        if destination.size != source.size:
            raise ValueError('Partitioned copy sizes differ')
        destination.bytes()[:] = source.bytes()
        return Completed()
    return cl.enqueue_copy(queue, destination, source)


class PartitionedPool:
    partitioned = True

    def __init__(self, queues, weights, min_cells, memory_fraction=0.8):
        self.queues = queues
        self.weights = weights
        self.min_cells = min_cells
        self.stats = {}
        self.memory_stats = {}
        self.memory_fraction = memory_fraction
        self.cell_order = None

    def repartition(self, centers):
        """Recursive spatial bisection balances population, not cell ID ranges."""
        xyz = np.asarray(centers).reshape(-1).view(np.float32).reshape(-1, 4)[:, :3]
        self.cell_order = spatial_order(xyz, self.weights)

    def ownership(self, count):
        order = self.cell_order if self.cell_order is not None and len(self.cell_order) == count else np.arange(count, dtype=np.int32)
        # Keep all selected GPUs active even below the compute-only threshold:
        # collapsing to GPU 0 could reintroduce the memory capacity limit.
        return [np.sort(order[start:start+size]) for start, size in partitions(count, self.weights)]

    def record(self, name, owned, trailing=1):
        counts = self.stats.setdefault(name, [0]*len(self.queues))
        for index, ids in enumerate(owned):
            counts[index] += len(ids)*trailing

    def check_memory(self, device, sizes):
        hardware = self.queues[device].device
        largest = max(sizes, default=0)
        if largest > hardware.max_mem_alloc_size or sum(sizes) > hardware.global_mem_size*self.memory_fraction:
            raise MemoryError('Partition on GPU %s needs %d bytes (largest buffer %d); '
                              'reduce cells/contact capacity or add GPUs. The spatial halo may be too large.' %
                              (hardware.name, sum(sizes), largest))

    def clear_cache(self):
        # No persistent device replicas; completed stages release scratch buffers.
        for queue in self.queues:
            queue.finish()


def spatial_order(xyz, weights):
    """Balanced recursive bisection; preserves global cell IDs through migration."""
    ids = np.arange(len(xyz), dtype=np.int32)
    sizes = [size for _, size in partitions(len(ids), weights)]
    def split(indices, counts):
        if len(counts) == 1 or not len(indices):
            return indices
        middle = len(counts)//2
        axis = int(np.argmax(np.ptp(xyz[indices], axis=0)))
        # Stable tie-breaking makes equal-position cells deterministic.
        indices = indices[np.argsort(xyz[indices, axis], kind='stable')]
        boundary = sum(counts[:middle])
        return np.concatenate((split(indices[:boundary], counts[:middle]),
                               split(indices[boundary:], counts[middle:])))
    return split(ids, sizes)


@dataclass
class Region:
    ids: object
    width: int = 1  # typed elements per global row

    def __post_init__(self):
        indices = np.unique(np.asarray(self.ids, dtype=np.int64))
        if len(indices) and (indices[0] < 0 or indices[-1] > np.iinfo(np.int32).max):
            raise ValueError('Partition indices exceed the OpenCL int index range')
        self.ids = indices.astype(np.int32)
        if self.width < 1 or np.any(self.ids < 0):
            raise ValueError('Invalid partition region')

    def pack(self, buffer, itemsize):
        rowsize = self.width*itemsize
        raw = buffer.bytes()
        if len(raw) % rowsize:
            # Legacy scratch buffers can be reinterpreted (float <-> float8).
            # Ignore only unused trailing storage, never an accessed row.
            raw = raw[:len(raw)//rowsize*rowsize]
        rows = raw.reshape(-1, rowsize)
        if len(self.ids) and self.ids[-1] >= len(rows):
            raise ValueError('Partition region exceeds host buffer capacity')
        return np.ascontiguousarray(rows[self.ids])

    def commit(self, buffer, itemsize, packed, owned):
        owned = np.asarray(owned, dtype=np.int32)
        positions = np.searchsorted(self.ids, owned)
        if len(owned) and (np.any(positions >= len(self.ids)) or np.any(self.ids[positions] != owned)):
            raise ValueError('Output ownership is not contained in its region')
        raw = buffer.bytes()
        width = self.width*itemsize
        raw[:len(raw)//width*width].reshape(-1, width)[owned] = packed[positions]


@dataclass
class Argument:
    name: str
    dtype: str
    pointer: bool

    @property
    def itemsize(self):
        return {'int': 4, 'float': 4, 'float4': 16, 'float8': 32}[self.dtype]


def kernel_parts(source):
    """Extract the repository's simple kernel declarations (not a C parser)."""
    result = {}
    for match in re.finditer(r'__kernel\s+void\s+(\w+)\s*\((.*?)\)\s*\{', source, re.S):
        depth, end = 1, match.end()
        while depth and end < len(source):
            depth += (source[end] == '{') - (source[end] == '}')
            end += 1
        if depth:
            raise ValueError('Unclosed OpenCL kernel')
        args = []
        for declaration in match.group(2).split(','):
            parsed = re.fullmatch(r'\s*(?:__global\s+)?(?:const\s+)?(\w+)\s*(\*?)\s*(\w+)\s*', declaration)
            if not parsed:
                raise ValueError('Unsupported partitioned kernel declaration: '+declaration)
            dtype, pointer, name = parsed.groups()
            args.append(Argument(name, dtype, bool(pointer)))
        result[match.group(1)] = (match.start(), end, match.group(2), source[match.end():end-1], args)
    return result


LOOKUP_SOURCE = r'''
// Missing halo entries abort the stage; never silently read an unrelated cell.
long cm_lookup(__global const int *ids, const int count, const int width,
               const long index, __global unsigned int *error) {
    if (index < 0) { atomic_or(error, 1u); return 0; }
    long row = index / width;
    int lo=0, hi=count;
    while (lo<hi) { int mid=lo+(hi-lo)/2; if (ids[mid]<row) lo=mid+1; else hi=mid; }
    if (lo>=count || ids[lo]!=row) { atomic_or(error, 1u); return 0; }
    return (long)lo*width + index%width;
}
'''


def rewrite_accesses(body, names):
    """Remap nested subscripts, plus the cell-local pointer offsets in integrators."""
    pattern = re.compile(r'\b('+'|'.join(map(re.escape, names))+r')\s*\[')
    def lookup(name, index):
        return 'cm_lookup(cm_%s_ids, cm_%s_count, cm_%s_width, (long)(%s), cm_error)' % (name, name, name, index)
    def rewrite(text):
        pieces, offset = [], 0
        while True:
            match = pattern.search(text, offset)
            if not match:
                pieces.append(text[offset:])
                break
            pieces.append(text[offset:match.start()])
            end, depth = match.end(), 1
            while end < len(text) and depth:
                depth += (text[end] == '[') - (text[end] == ']')
                end += 1
            if depth:
                raise ValueError('Unclosed array subscript')
            index = rewrite(text[match.end():end-1])
            pieces.append(match.group(1)+'['+lookup(match.group(1), index)+']')
            offset = end
        return ''.join(pieces)
    body = rewrite(body)
    for name in names:
        # Repository kernels pass a pointer to the start of one cell's block to
        # user rate functions. Those functions keep their ordinary local indices.
        body = re.sub(r'\b'+re.escape(name)+r'\s*\+\s*([A-Za-z_]\w*)',
                      lambda m: name+' + '+lookup(name, m.group(1)), body)
    return body


def partition_source(source):
    pieces, offset = [LOOKUP_SOURCE], 0
    for name, (start, end, declarations, body, args) in kernel_parts(source).items():
        pieces.append(source[offset:start])
        names = [arg.name for arg in args if arg.pointer]
        # Strip comments before mechanical access rewriting. Do not rewrite
        # user functions; only the audited built-in kernels are partitioned.
        body = re.sub(r'/\*.*?\*/|//[^\n]*', '', body, flags=re.S)
        if name.startswith('find_') and name.endswith('_contacts'):
            body = body.replace('k++;', 'if (k >= max_contacts) { atomic_or(cm_error, 2u); return; } k++;')
        if name == 'collect_tos':
            body = body.replace('cell_tos[i*max_contacts+k] =',
                                'if (k >= max_contacts) { atomic_or(cm_error, 2u); return; } cell_tos[i*max_contacts+k] =')
        body = rewrite_accesses(body, names)
        body = body.replace('get_global_id(0)', 'cm_owned[get_global_id(0)]')
        extra = ['__global const int *cm_owned', '__global unsigned int *cm_error']
        for pointer in names:
            extra.extend(['__global const int *cm_%s_ids' % pointer,
                          'int cm_%s_count' % pointer, 'int cm_%s_width' % pointer])
        pieces.append('__kernel void %s(%s, %s) {\n%s\n}' % (name, declarations, ', '.join(extra), body))
        offset = end
    pieces.append(source[offset:])
    return ''.join(pieces)


def _integers(buffer):
    return buffer.array.reshape(-1).view(np.int32).reshape(-1)


def _union(*arrays):
    return np.unique(np.concatenate([np.asarray(a, dtype=np.int64).reshape(-1) for a in arrays])).astype(np.int32)


def _contact_indices(cells, counts, stride):
    lengths = counts[cells]
    if np.any(lengths < 0) or np.any(lengths > stride):
        raise ValueError('Contact count exceeds max_contacts')
    return (cells[:, None]*stride+np.arange(stride)[None, :])[np.arange(stride)[None, :] < lengths[:, None]]


def neighbor_regions(values, owned):
    """Exact grid bins visited by the existing contact kernels, including ends."""
    sx = int(values['grid_x_max'])-int(values['grid_x_min'])
    sy = int(values['grid_y_max'])-int(values['grid_y_min'])
    nsquares, ncells = int(values['n_sqs']), int(values['n_cells'])
    sqs = _integers(values['sqs'])[owned]
    if sx <= 0 or sy <= 0 or nsquares != sx*sy or np.any(sqs < 0) or np.any(sqs >= nsquares):
        raise ValueError('Invalid spatial grid for halo exchange')
    bins = []
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            y, x = sqs//sx+dy, sqs%sx+dx
            valid = (y >= 0) & (y < sy) & (x >= 0) & (x < sx)
            bins.append(y[valid]*sx+x[valid])
    bins = _union(*bins)
    starts = _integers(values['sq_inds'])
    ranks = [np.arange(starts[b], starts[b+1] if b < nsquares-1 else ncells, dtype=np.int32) for b in bins]
    ranks = _union(*ranks) if ranks else np.empty(0, np.int32)
    neighbors = np.unique(_integers(values['sorted_ids'])[ranks])
    if np.any(neighbors < 0) or np.any(neighbors >= ncells):
        raise ValueError('Invalid cell ID in spatial index')
    boundaries = _union(bins, bins[bins < nsquares-1]+1)
    return neighbors, ranks, boundaries


def plan_regions(name, arguments, args, owned, outputs, family):
    """Audited per-buffer read sets; no full-buffer fallback for unknown stages."""
    values = {arg.name: value for arg, value in zip(arguments, args)}
    regions = {arg.name: Region(owned) for arg in arguments if arg.pointer}
    # Every output's row stride comes from the same audited layout as the
    # compute-only backend, expressed in the actual kernel pointer type.
    for index, width, offset in outputs:
        if offset or width % arguments[index].itemsize:
            raise ValueError('Unsupported partition output layout')
        regions[arguments[index].name] = Region(owned, max(1, width//arguments[index].itemsize))

    def assign(names, ids, width=1):
        for key in names.split():
            if key in regions:
                regions[key] = Region(ids, width)

    if family == 'physics':
        if name not in ('bin_cells', 'find_plane_contacts', 'find_sphere_contacts', 'find_contacts',
                        'collect_tos', 'build_matrix', 'calculate_Bx', 'calculate_BTBx',
                        'calculate_Mx', 'calculate_Minv_x', 'predict', 'integrate', 'add_impulse'):
            raise ValueError('Unaudited physics kernel '+name)
        stride = int(values.get('max_contacts', 1))
        assign('frs tos dists pts norms reldists stiff overlap fr_ents to_ents Bx cell_tos', owned, stride)
        # calculate_Mx has a float8 output despite its shared scratch buffer
        # being allocated as float in the legacy model. Typed byte packing handles it.
        if name in ('find_plane_contacts', 'find_sphere_contacts'):
            prefix = 'plane' if name == 'find_plane_contacts' else 'sphere'
            ids = np.arange(int(values['n_'+prefix+'s']), dtype=np.int32)
            assign(prefix+'_pts '+prefix+'_norms '+prefix+'_coeffs '+prefix+'_rads', ids)
        elif name in ('find_contacts', 'collect_tos'):
            neighbors, ranks, boundaries = neighbor_regions(values, owned)
            assign('sorted_ids', ranks)
            assign('sq_inds', boundaries)
            if name == 'find_contacts':
                assign('centers dirs lens rads', _union(owned, neighbors))
            else:
                assign('n_cts', neighbors)
                # The reverse-contact scan reads only occupied target slots.
                # frs is an unused legacy argument in this kernel.
                assign('tos', _contact_indices(neighbors, _integers(values['n_cts']), stride))
                assign('frs', [])
        elif name in ('build_matrix', 'calculate_Bx'):
            frs, tos = _integers(values['frs']), _integers(values['tos'])
            if name == 'build_matrix':
                indices = _contact_indices(owned, _integers(values['n_cts']), stride)
            else:
                indices = (owned[:, None]*stride+np.arange(stride)[None, :]).reshape(-1)
                indices = indices[(frs[indices] != 0) | (tos[indices] != 0)]
            endpoints = _union(frs[indices], tos[indices][tos[indices] >= 0])
            assign('centers dirs lens rads deltap', endpoints)
        elif name == 'calculate_BTBx':
            outgoing = _contact_indices(owned, _integers(values['n_cts']), stride)
            incoming_slots = _contact_indices(owned, _integers(values['n_cell_tos']), stride)
            incoming = _integers(values['cell_tos'])[incoming_slots]
            incoming = incoming[incoming >= 0]
            # Exchange only the referenced contact entries, not the complete
            # contact table of neighboring cells, during transpose products.
            assign('fr_ents', outgoing)
            assign('to_ents', incoming)
            assign('Bx', _union(outgoing, incoming))
    elif family in ('euler', 'signal'):
        nspecies, nsignals = int(values.get('numSpecies', 1)), int(values.get('numSignals', 1))
        implicit = 'Implicit' in name
        assign('cellSpecLevels specRate specRates', owned, nspecies)
        assign('cellSignalLevels', owned, nsignals*(2 if implicit else 1))
        assign('weights indices', owned, 8)
        assign('sigRates', owned, 8*nsignals)
        if name.startswith('setCellSig'):
            assign('levels', owned, nsignals*(2 if implicit else 1))
            indices = _integers(values['indices']).reshape(-1, 8)[owned].reshape(-1)
            grid_size = int(values['gridTotalSize'])
            if np.any(indices < 0) or np.any(indices >= grid_size):
                raise ValueError('Signal-grid halo contains invalid indices')
            grid_ids = np.unique((indices[:, None]+np.arange(nsignals)[None, :]*grid_size).reshape(-1))
            assign('grid transport', grid_ids)
    elif family not in ('elementwise', 'fixed'):
        raise ValueError('Unknown partitioned kernel family '+family)
    if family == 'physics' and 'max_contacts' in values:
        for count_name in ('n_cts', 'n_cell_tos'):
            if count_name in regions:
                counts = _integers(values[count_name])[regions[count_name].ids]
                if np.any(counts < 0) or np.any(counts > int(values['max_contacts'])):
                    raise ValueError('Contact count exceeds max_contacts')
    return [regions.get(arg.name) for arg in arguments]


class PartitionedProgram:
    def __init__(self, source, pool, layouts, label='', family=None):
        self.pool, self.layouts, self.label = pool, layouts, label
        self.family = family or ('physics' if label == 'physics.' else 'fixed' if label == 'fixed.'
                                 else 'euler' if label == 'CLEulerIntegrator.' else 'signal')
        self.parts = kernel_parts(source)
        if set(self.parts) != set(layouts):
            raise ValueError('Partitioned execution requires a layout for every kernel')
        transformed = partition_source(source)
        self.programs = [cl.Program(queue.context, transformed).build(cache_dir=False) for queue in pool.queues]
        self.kernels = {}

    def __getattr__(self, name):
        if name not in self.layouts:
            raise AttributeError('No partitioned layout for '+name)
        if name in self.kernels:
            return self.kernels[name]
        kernels = [getattr(program, name) for program in self.programs]
        arguments = self.parts[name][4]
        def run(queue, shape, local, *args, **kwargs):
            if queue != self.pool.queues[0] or local is not None or kwargs:
                raise ValueError('Partitioned kernels require a plain launch on the simulator queue')
            if any(n == 0 for n in shape):
                return Completed()
            if len(args) != len(arguments) or any(arg.pointer and not isinstance(value, HostBuffer)
                                                  for arg, value in zip(arguments, args)):
                raise TypeError('Partitioned kernels require partitioned array handles')
            owned = self.pool.ownership(shape[0])
            outputs = self.layouts[name](args)
            prepared = []
            usage = []
            for device, ids in enumerate(owned):
                if not len(ids):
                    usage.append({'bytes': 0, 'owned': 0, 'buffers': {}})
                    continue
                regions = plan_regions(name, arguments, args, ids, outputs, self.family)
                packed = [region.pack(value, arg.itemsize) if region is not None else None
                          for arg, value, region in zip(arguments, args, regions)]
                sizes = [ids.nbytes, 4]
                buffers = {}
                for arg, region, data in zip(arguments, regions, packed):
                    if region is not None:
                        sizes.extend([max(region.width*arg.itemsize, data.nbytes), max(4, region.ids.nbytes)])
                        buffers[arg.name] = {'rows': len(region.ids), 'bytes': data.nbytes,
                                             'row_width': region.width}
                self.pool.check_memory(device, sizes)
                usage.append({'bytes': sum(sizes), 'owned': len(ids), 'buffers': buffers})
                prepared.append((device, ids, regions, packed))
            # Preflight every device before allocating or launching anything.
            workers = []
            try:
                for device, ids, regions, packed in prepared:
                    worker_queue = self.pool.queues[device]
                    raw, metadata, references = [], [], []
                    def upload(array):
                        value = np.ascontiguousarray(array)
                        if not value.nbytes:
                            value = np.zeros(1, np.int32)
                        buffer = cl.Buffer(worker_queue.context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                                           hostbuf=value)
                        references.append(buffer)
                        return buffer
                    owner_dev, error_dev = upload(ids), upload(np.zeros(1, np.uint32))
                    for arg, value, region, data in zip(arguments, args, regions, packed):
                        if region is None:
                            raw.append(value)
                        else:
                            # A failed lookup returns index zero while setting
                            # the stage error flag. Even an empty halo needs a
                            # full typed row so that this error path is in bounds.
                            device_data = data if data.nbytes else np.zeros((1, region.width*arg.itemsize), np.uint8)
                            raw.append(upload(device_data))
                            metadata.extend([upload(region.ids), np.int32(len(region.ids)), np.int32(region.width)])
                    event = kernels[device](worker_queue, (len(ids),)+tuple(shape[1:]), None,
                                   *raw, owner_dev, error_dev, *metadata)
                    worker_queue.flush()
                    workers.append((device, ids, regions, packed, raw, error_dev, event, references))
                for device, ids, regions, packed, raw, error_dev, event, references in workers:
                    event.wait()
                    error = np.zeros(1, np.uint32)
                    cl.enqueue_copy(self.pool.queues[device], error, error_dev).wait()
                    if error[0]:
                        raise RuntimeError('Partitioned kernel %s failed on GPU %d: %s' %
                                           (name, device, 'contact capacity exceeded' if error[0] & 2 else 'missing halo index'))
                    for index, width, _ in outputs:
                        if width:
                            cl.enqueue_copy(self.pool.queues[device], packed[index], raw[index]).wait()
                # Commit only after all workers/halos succeed, keeping host
                # state unchanged if a launch reports an error.
                for device, ids, regions, packed, raw, error_dev, event, references in workers:
                    for index, width, _ in outputs:
                        if width:
                            regions[index].commit(args[index], arguments[index].itemsize, packed[index], ids)
            finally:
                for worker_queue in self.pool.queues:
                    worker_queue.finish()
            self.pool.record(self.label+name, owned, math.prod(shape[1:]))
            self.pool.memory_stats[self.label+name] = usage
            return Completed()
        self.kernels[name] = run
        return run


class PartitionedElementwise:
    def __init__(self, pool, arguments, operation, name):
        # These built-in vector kernels use row-local pointer subscripts only.
        declarations = []
        for declaration in arguments.split(','):
            declarations.append(('__global ' if '*' in declaration else '')+declaration.strip())
        source = '__kernel void %s(%s) { int i=get_global_id(0); %s; }' % (name, ','.join(declarations), operation)
        self.args = kernel_parts(source)[name][4]
        self.program = PartitionedProgram(source, pool, {name: lambda a: [(0, self.args[0].itemsize, 0)]},
                                          family='elementwise')
        self.pool, self.name = pool, name

    def __call__(self, *args, **kwargs):
        if kwargs:
            raise ValueError('Partitioned elementwise kernels expect positional arguments')
        raw = [arg.data if spec.pointer else {'float': np.float32, 'int': np.int32}[spec.dtype](arg)
               for arg, spec in zip(args, self.args)]
        return getattr(self.program, self.name)(self.pool.queues[0], (args[0].size,), None, *raw)


class PartitionedReduction:
    def __init__(self, original, pool, dtype, name, factory=None):
        self.pool, self.dtype, self.name = pool, dtype, name
        self.originals = [original if queue.context == pool.queues[0].context else factory(queue.context)
                          for queue in pool.queues]

    def __call__(self, *args, **kwargs):
        import pyopencl.array as arrays
        if kwargs or any(not isinstance(a, PartitionedArray) or not a.flags.c_contiguous or a.size != args[0].size for a in args):
            raise ValueError('Partitioned reduction requires equally-sized contiguous arrays')
        owned = self.pool.ownership(args[0].size)
        local = [[np.ascontiguousarray(arg.array.reshape(-1)[ids]) for arg in args] for ids in owned]
        usage = []
        for device, inputs in enumerate(local):
            sizes = [a.nbytes for a in inputs]
            # The generated tree reduction uses temporary partial sums; reserve
            # an extra input-sized amount for a conservative memory preflight.
            if len(owned[device]):
                self.pool.check_memory(device, sizes+[max(sizes)])
            usage.append({'bytes': sum(sizes)+max(sizes, default=0), 'owned': len(owned[device])})
        workers = []
        try:
            for device, (ids, inputs) in enumerate(zip(owned, local)):
                if len(ids):
                    queue = self.pool.queues[device]
                    device_arrays = [arrays.to_device(queue, value) for value in inputs]
                    result = self.originals[device](*device_arrays, queue=queue)
                    queue.flush()
                    workers.append((result, device_arrays))
            values = [result.get().item() for result, _ in workers]
        finally:
            for queue in self.pool.queues:
                queue.finish()
        total = math.fsum(values) if np.dtype(self.dtype).kind == 'f' else sum(values)
        self.pool.record(self.name, owned)
        self.pool.memory_stats[self.name] = usage
        return PartitionedArray(np.asarray(total, dtype=self.dtype), self.pool.queues[0])
