"""Partition-resident CG: static matrix and vectors stay on device per solve.

Only endpoint-vector halos and reduction scalars cross the host during CG.
Separate OpenCL contexts require host-staged halo exchange (not peer-to-peer).
The final displacement is committed only after successful finite readback.
"""
import math
from pathlib import Path

import numpy as np


SOURCE = r'''
float rd_dot(float8 a, float8 b) {
    return dot(a.s0123,b.s0123)+dot(a.s4567,b.s4567);
}
__kernel void rd_bx(__global const int *fr, __global const int *to,
    __global const int *valid, __global const float8 *f, __global const float8 *t,
    __global const float8 *p, __global float *bx) {
    int i=get_global_id(0), b=to[i];
    float8 pa=p[fr[i]], pb=b<0?(float8)(0):p[b];
    float8 ti=b<0?(float8)(0):t[i];
    bx[i]=valid[i] ? (dot(f[i].s0123,pa.s0123)-dot(ti.s0123,pb.s0123))
                    +(dot(f[i].s4567,pa.s4567)-dot(ti.s4567,pb.s4567)) : 0.f;
}
__kernel void rd_transpose(__global const int *offsets, __global const int *contacts,
    __global const int *signs, __global const float8 *f, __global const float8 *t,
    __global const float *bx, __global float8 *result) {
    int i=get_global_id(0); float8 value=0.f;
    for(int k=offsets[i]; k<offsets[i+1]; k++) {
        int c=contacts[k];
        if(signs[k]>0) value+=f[c]*bx[c]; else value-=t[c]*bx[c];
    }
    result[i]=value;
}
__kernel void rd_regularize(float scale, __global float8 *ap, __global const float8 *mx) {
    int i=get_global_id(0); ap[i]+=scale*mx[i];
}
__kernel void rd_update(float alpha, __global float8 *x, __global float8 *r,
    __global const float8 *p, __global const float8 *ap) {
    int i=get_global_id(0); x[i]+=alpha*p[i]; r[i]-=alpha*ap[i];
}
__kernel void rd_direction(float beta, __global float8 *p, __global const float8 *r) {
    int i=get_global_id(0); p[i]=r[i]+beta*p[i];
}
__kernel void rd_gather(__global const int *indices, __global const float8 *p,
    __global float8 *packed) { int i=get_global_id(0); packed[i]=p[indices[i]]; }
'''


def plan_partitions(owned, n, stride, counts, incoming_counts, incoming, fr, to):
    """Construct compact contact rows and explicit cross-device gather routes."""
    order = np.concatenate(owned)
    if len(order) != n or not np.array_equal(np.sort(order), np.arange(n)):
        raise ValueError('Resident solver ownership must cover each cell exactly once')
    for values in (counts[:n], incoming_counts[:n]):
        if np.any(values < 0) or np.any(values > stride):
            raise ValueError('Invalid resident contact count')
    packets = []
    for ids in owned:
        row_contacts, signs, offsets = [], [], [0]
        for cell in ids:
            outgoing = list(range(int(cell)*stride, int(cell)*stride+int(counts[cell])))
            backward = incoming[cell*stride:cell*stride+incoming_counts[cell]]
            backward = backward[backward >= 0].tolist()
            row_contacts.extend(outgoing + backward)
            signs.extend([1]*len(outgoing)+[-1]*len(backward))
            offsets.append(len(row_contacts))
        contacts = np.unique(np.asarray(row_contacts, dtype=np.int32))
        if len(contacts) and (contacts[0] < 0 or contacts[-1] >= n*stride):
            raise ValueError('Resident contact reference outside active population')
        endpoints = np.unique(np.concatenate((fr[contacts], to[contacts][to[contacts] >= 0])))
        if len(endpoints) and (endpoints[0] < 0 or endpoints[-1] >= n):
            raise ValueError('Resident contact endpoint outside active population')
        halo = np.setdiff1d(endpoints, ids)
        # Group remote cells by source device so every received block is contiguous.
        remote = [np.intersect1d(halo, source) for source in owned]
        columns = np.concatenate([ids]+remote)
        lookup = np.full(n, -1, np.int32)
        lookup[columns] = np.arange(len(columns), dtype=np.int32)
        packets.append(dict(owned=ids, columns=columns, remote=remote,
                            contact_ids=contacts, offsets=np.asarray(offsets, np.int32),
                            contacts=np.searchsorted(contacts, row_contacts).astype(np.int32),
                            signs=np.asarray(signs, np.int32), fr=lookup[fr[contacts]],
                            to=np.where(to[contacts] < 0, -1, lookup[np.maximum(to[contacts], 0)]).astype(np.int32),
                            valid=((fr[contacts] != 0) | (to[contacts] != 0)).astype(np.int32)))
    return packets


class ResidentCG:
    def __init__(self, pool):
        import pyopencl as cl
        from pyopencl.reduction import ReductionKernel
        self.pool = pool
        # Use the original inertia implementation and calculate_Mx kernel.
        physics = (Path(__file__).parent / 'Biophysics/BacterialModels/CLBacterium.cl').read_text()
        self.programs = [cl.Program(q.context, physics+'\n'+SOURCE).build() for q in pool.queues]
        self.dots = [ReductionKernel(q.context, np.float32, neutral='0', reduce_expr='a+b',
                      map_expr='dot(x[i].s0123,y[i].s0123)+dot(x[i].s4567,y[i].s4567)',
                      arguments='__global float8 *x, __global float8 *y') for q in pool.queues]
        self.stats = {}

    def solve(self, phys):
        import pyopencl as cl
        import pyopencl.array as arrays
        from pyopencl.array import vec
        pool, n, stride = self.pool, phys.n_cells, phys.max_contacts
        if not n:
            return (0, 0.0)
        if not math.isfinite(phys.gamma) or phys.gamma <= 0 or phys.cgs_tol <= 0:
            raise ValueError('Resident CG requires positive gamma and solver tolerance')
        def host(name, dtype=np.float32, width=1):
            return getattr(phys, name).array.view(dtype).reshape(-1, width) if width > 1 else getattr(phys, name).array.view(dtype).reshape(-1)
        owned = pool.ownership(n)
        plans = plan_partitions(owned, n, stride, host('cell_n_cts_dev', np.int32),
                                host('n_cell_tos_dev', np.int32), host('cell_tos_dev', np.int32),
                                host('ct_frs_dev', np.int32), host('ct_tos_dev', np.int32))
        stats = dict(iterations=0, halo_bytes=0, scalar_bytes=0, solution_bytes=0,
                     owned_cells=[len(ids) for ids in owned], planned_bytes=[])
        self.stats.clear()
        self.stats.update(stats)
        packets = []
        for plan in plans:
            ct, ids = plan['contact_ids'], plan['owned']
            data = {key: plan[key] for key in ('offsets', 'contacts', 'signs', 'fr', 'to', 'valid')}
            data.update(f=host('fr_ents_dev', width=8)[ct], t=host('to_ents_dev', width=8)[ct],
                        rhs_contact=host('ct_reldists_dev')[ct], dirs=host('cell_dirs_dev', width=4)[ids],
                        lens=host('cell_lens_dev')[ids], rads=host('cell_rads_dev')[ids])
            packets.append(data)
        routes = []
        for destination, plan in enumerate(plans):
            offset = len(plan['owned'])
            for source, ids in enumerate(plan['remote']):
                if len(ids):
                    routes.append(dict(source=source, destination=destination, offset=offset,
                                       indices=np.searchsorted(owned[source], ids).astype(np.int32)))
                    offset += len(ids)
        # Account for resident vectors, halo staging, static matrix and reductions
        # before any allocation. Drop stage caches while a resident solve owns VRAM.
        pool.clear_cache()
        for device, (plan, packet) in enumerate(zip(plans, packets)):
            size, columns, contacts = len(plan['owned']), len(plan['columns']), len(plan['contact_ids'])
            sizes = [max(4, value.nbytes) for value in packet.values()]
            sizes += [max(32, columns*32)] + [max(32, size*32)]*5 + [max(4, contacts*4)]
            for route in routes:
                if route['source'] == device:
                    sizes += [route['indices'].nbytes, len(route['indices'])*32]
            sizes += [max(32, size*32)]  # conservative reduction temporary reserve
            if size:
                pool.check_memory(device, sizes)
            self.stats['planned_bytes'].append(sum(sizes) if size else 0)
        workers = []
        try:
            for device, (plan, packet) in enumerate(zip(plans, packets)):
                q, size = pool.queues[device], len(plan['owned'])
                if not size:
                    workers.append(None)
                    continue
                buffers = {key: arrays.to_device(q, np.ascontiguousarray(value) if value.size else np.zeros(1, np.int32))
                           for key, value in packet.items()}
                for key in ('x', 'r', 'ap', 'mx', 'rhs'):
                    buffers[key] = arrays.zeros(q, size, vec.float8)
                buffers['p'] = arrays.zeros(q, len(plan['columns']), vec.float8)
                buffers['bx'] = arrays.zeros(q, max(1, len(plan['contact_ids'])), np.float32)
                workers.append(buffers)
            for route in routes:
                q = pool.queues[route['source']]
                route['gpu_indices'] = arrays.to_device(q, route['indices'])
                route['gpu_values'] = arrays.empty(q, len(route['indices']), vec.float8)
                route['host_values'] = np.empty(len(route['indices']), vec.float8)

            def transpose(device, source, result):
                b, q, program = workers[device], pool.queues[device], self.programs[device]
                program.rd_transpose(q, (len(owned[device]),), None,
                                     *[b[k].data for k in ('offsets', 'contacts', 'signs', 'f', 't')],
                                     b[source].data, b[result].data)

            def exchange():
                events = []
                for route in routes:
                    source = route['source']
                    q = pool.queues[source]
                    self.programs[source].rd_gather(q, (len(route['indices']),), None,
                        route['gpu_indices'].data, workers[source]['p'].data, route['gpu_values'].data)
                    events.append(cl.enqueue_copy(q, route['host_values'], route['gpu_values'].data, is_blocking=False))
                    q.flush()
                for route, event in zip(routes, events):
                    event.wait()
                    cl.enqueue_copy(pool.queues[route['destination']], workers[route['destination']]['p'].data,
                                    route['host_values'], device_offset=route['offset']*32).wait()
                    self.stats['halo_bytes'] += 2*route['host_values'].nbytes

            def dot(left, right):
                partials = []
                for device, b in enumerate(workers):
                    if b is not None:
                        q = pool.queues[device]
                        partials.append(self.dots[device](b[left][:len(owned[device])],
                                                         b[right][:len(owned[device])], queue=q))
                        q.flush()
                total = np.float32(math.fsum(float(value.get()) for value in partials))
                self.stats['scalar_bytes'] += 4*len(partials)
                pool.record('physics.dot', owned)
                if not np.isfinite(total):
                    raise FloatingPointError('Nonfinite resident CG reduction')
                return total

            for device, b in enumerate(workers):
                if b is not None:
                    transpose(device, 'rhs_contact', 'rhs')
                    q = pool.queues[device]
                    cl.enqueue_copy(q, b['r'].data, b['rhs'].data)
                    cl.enqueue_copy(q, b['p'].data, b['rhs'].data, byte_count=len(owned[device])*32)
            rsold = dot('r', 'r')
            residual = math.sqrt(float(rsold)/n)
            iterations = 0
            while residual >= phys.cgs_tol and iterations < n*7:
                exchange()
                for device, b in enumerate(workers):
                    if b is None:
                        continue
                    q, program = pool.queues[device], self.programs[device]
                    contacts = len(plans[device]['contact_ids'])
                    if contacts:
                        program.rd_bx(q, (contacts,), None,
                                      *[b[k].data for k in ('fr', 'to', 'valid', 'f', 't', 'p', 'bx')])
                    transpose(device, 'bx', 'ap')
                    program.calculate_Mx(q, (len(owned[device]),), None, np.float32(phys.muA),
                                         np.float32(phys.gamma), *[b[k].data for k in ('dirs', 'lens', 'rads', 'p', 'mx')])
                    program.rd_regularize(q, (len(owned[device]),), None, np.float32(1/phys.gamma), b['ap'].data, b['mx'].data)
                pool.record('physics.calculate_Bx', owned, stride)
                pool.record('physics.calculate_BTBx', owned)
                p_ap = dot('p', 'ap')
                if p_ap <= 0:
                    raise FloatingPointError('Resident CG breakdown: nonpositive p^T A p')
                alpha = np.float32(rsold/p_ap)
                for device, b in enumerate(workers):
                    if b is not None:
                        self.programs[device].rd_update(pool.queues[device], (len(owned[device]),), None,
                                                        alpha, *[b[k].data for k in ('x', 'r', 'p', 'ap')])
                pool.record('vecaddkx', owned)
                rsnew = dot('r', 'r')
                residual = math.sqrt(float(rsnew)/n)
                iterations += 1
                if residual < phys.cgs_tol:
                    break
                beta = np.float32(rsnew/rsold)
                for device, b in enumerate(workers):
                    if b is not None:
                        self.programs[device].rd_direction(pool.queues[device], (len(owned[device]),), None,
                                                           beta, b['p'].data, b['r'].data)
                rsold = rsnew
            # Stage the complete solution before committing any partition.
            results = [(ids, b['x'].get()) for ids, b in zip(owned, workers) if b is not None]
            if any(not np.isfinite(value.view(np.float32)).all() for _, value in results):
                raise FloatingPointError('Nonfinite resident CG solution')
            for ids, value in results:
                phys.deltap_dev.array[ids] = value
            self.stats.update(iterations=iterations, residual=residual, solution_bytes=n*32,
                              converged=residual < phys.cgs_tol)
            return iterations, residual
        finally:
            for q in pool.queues:
                q.finish()
