import argparse
import time 
import os
import sys

DEFAULT_OUT_ROOT = '~/cellmodeller'

def parse_args():
    '''
    Parse command line arguments for running a CellModeller model headlessly (without a GUI)on a server.
    '''

    parser = argparse.ArgumentParser(description='CellModeller headless run for server simulations')
    parser.add_argument("model", nargs='?', 
                        help="Path to model file fx. Project/anti_square_slow_server.py")
    parser.add_argument("--max-cells", type=int, default=10*int(1e4),
                        help="Maximum number of cells to simulate (default: 100,000)")
    # Change default when running on server
    parser.add_argument("--max-time", type=int, default=60,
                        help="Maximum simulation time in seconds (default: 60s i.e. 1 min)")
    parser.add_argument("--out-root", type=str, default=DEFAULT_OUT_ROOT,
                        help="Root directory for simulation output (default: %s)" % DEFAULT_OUT_ROOT)
    parser.add_argument("--pickle-steps", type=int, default=100,
                        help="Number of steps between saving simulation state (default: 100)")
    parser.add_argument("--dt", type=float, default=0.05,
                        help="Time step for simulation (default: 0.05)")
    parser.add_argument("--name", type=str, default=None,
                        help="Name of the output folder (default: <model>-<date>)")
    parser.add_argument("--platform", type=int, default=0,
                        help="OpenCL platform number (default: 0)")
    parser.add_argument("--device", type=int, default=0,
                        help="OpenCL device number (default: 0)")
    parser.add_argument("--log-every", type=int, default=50,
                        help="Print progress every N steps (default: 50)")
    parser.add_argument("--list-devices", action='store_true',
                        help="List available OpenCL platforms and devices and exit")
    return parser.parse_args()


def list_devices():
    '''
    List available OpenCL platforms and devices.
    '''
    import pyopencl as cl
    for p, platform in enumerate(cl.get_platforms()):
        print('Platform %d: %s' % (p, platform.name))
        for d, device in enumerate(platform.get_devices()):
            print('    Device %d: %s (%s)' % (d, device.name, cl.device_type.to_string(device.type)))

def main():
    args = parse_args()

    # Set the maximum number of cells as an environment variable for the model to read
    os.environ['MAX_CELLS'] = str(args.max_cells)

    # Early exit conditions
    if args.list_devices:
        list_devices()
        return
    if args.model is None:
        sys.exit('Error: a model file is required')

    # Set up the output directory
    out_root = os.path.abspath(os.path.expanduser(args.out_root))

    # Check if the output directory exists i.e. contains a 'data' folder, if not exit with an error message
    if not os.path.isdir(os.path.join(out_root, 'data')):
        sys.exit('Output folder %s does not exist - create it first.' % os.path.join(out_root, 'data'))
    os.environ['CMPATH'] = out_root

    # Import the model by module name, with its folder on sys.path
    path, filename = os.path.split(os.path.abspath(args.model))
    modname = os.path.splitext(filename)[0]
    sys.path.insert(0, path)


    # Initialize the simulator with the specified parameters
    from CellModeller.Simulator import Simulator
    sim = Simulator(modname, args.dt, outputDirName=args.name,
                    clPlatformNum=args.platform, clDeviceNum=args.device, saveOutput=True)
    
    # Set the number of steps between saving simulation state if specified (should this not be coordinated with the model file?)
    if args.pickle_steps is not None:
        sim.pickleSteps = args.pickle_steps

    # Stopping point: stay a margin below the preallocated max_cells. One step can add ~1-3% more
    # cells, and CLBacterium silently skips divisions once its buffers are full.
    cell_limit = args.max_cells - max(256, int(0.05 * args.max_cells))

    print("Starting simulation with max_cells=%d (stopping at %d cells), max_time=%d seconds"
          % (args.max_cells, cell_limit, args.max_time))
    print("Writing output to: %s" % sim.outputDirPath)


    # Run the simulation loop until the stopping conditions are met
    start = time.time()
    while True:
        n_cells = len(sim.cellStates)
        time_elapsed = time.time() - start

        # Check if the maximum amount of time has been reached
        if args.max_time is not None and time_elapsed >= args.max_time:
            print("Reached maximum simulation time of %d seconds." % args.max_time)
            break

        # Check if the maximum number of cells has been reached
        if n_cells >= cell_limit:
            print("Reached maximum number of cells: %d" % n_cells)
            break

        # Step the simulation forward
        sim.step()

        # Log progress every N steps as specified by the user
        if sim.stepNum % args.log_every == 0:
            print('step %6d   cells %7d   elapsed %8.1f s' % (sim.stepNum, len(sim.cellStates), time.time() - start), flush=True)   

    # Final log message indicating the simulation has completed
    print('Done: %d steps, %d cells, %.1f s' % (sim.stepNum, len(sim.cellStates), time.time() - start), flush=True)


if __name__ == '__main__':
    main()

