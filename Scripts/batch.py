import os
import sys
import subprocess
import string
import shutil

from CellModeller.Simulator import Simulator

max_cells = 50000
cell_buffer = 256

def simulate(modfilename, platform, device, steps=50, devices=None,
             device_weights=None, multi_gpu_min_cells=1024, gpu_memory="partitioned",
             stop_cells=max_cells):
    (path,name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, 0.025, clPlatformNum=platform, clDeviceNum=device, saveOutput=True,
                    clDeviceNums=devices, clDeviceWeights=device_weights,
                    clMultiGPUMinCells=multi_gpu_min_cells, clMultiGPUMemory=gpu_memory)
    if stop_cells <= cell_buffer:
        raise ValueError('Stop-cell limit must exceed the cell buffer (%d)' % cell_buffer)
    capacity = getattr(sim.phys, 'max_cells', stop_cells)
    if stop_cells > capacity:
        raise ValueError('The model max_cells (%d) must be at least --max-cells (%d)' % (capacity, stop_cells))
    while len(sim.cellStates) < stop_cells-cell_buffer:
        sim.step()

def main():
    import argparse
    import pyopencl as cl
    parser = argparse.ArgumentParser(description='Run a CellModeller simulation')
    parser.add_argument('model', help='Model Python file')
    parser.add_argument('platform', type=int, nargs='?')
    parser.add_argument('device', type=int, nargs='?', default=0)
    parser.add_argument('--devices', help="Comma-separated device indices, or 'auto' for all GPUs on the platform")
    parser.add_argument('--device-weights', help='Comma-separated positive workload weights in device order')
    parser.add_argument('--multi-gpu-min-cells', type=int, default=1024)
    parser.add_argument('--gpu-memory', choices=['partitioned', 'replicated'], default='partitioned')
    parser.add_argument('--max-cells', type=int, default=max_cells,
                        help='Stop population limit; increase max_cells in the model to match')
    args = parser.parse_args()
    try:
        devices = args.devices
        if devices is not None and devices != 'auto':
            devices = [int(value) for value in devices.split(',')]
        weights = None if args.device_weights is None else [float(value) for value in args.device_weights.split(',')]
    except ValueError:
        parser.error('Device indices must be integers and weights must be numbers')
    platnum = args.platform
    devnum = args.device
    if platnum is None:
        platforms = cl.get_platforms()
        for i, platform in enumerate(platforms):
            print('%d: %s' % (i, platform))
        platnum = int(input('Platform number: '))
        if devices is None:
            for i, device in enumerate(platforms[platnum].get_devices()):
                print('%d: %s' % (i, device))
            devnum = int(input('Device number: '))
    simulate(args.model, platnum, devnum, devices=devices,
             device_weights=weights, multi_gpu_min_cells=args.multi_gpu_min_cells,
             gpu_memory=args.gpu_memory, stop_cells=args.max_cells)

# Make sure we are running as a script
if __name__ == "__main__": 
    main()
