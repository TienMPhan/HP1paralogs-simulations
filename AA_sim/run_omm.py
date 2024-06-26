'''
  - Production run in OpenMM with files loaded from Gromacs.
  - The simulation will be restarted from the equilibration step via loadCheckpoint().
  - Checkpoint will be save for restart.
  - 2 places should be changed to run on different systems:
        i)  top = load_file('filename.top', xyz='filename.gro'
        ii) pro = load('filename.gro')
  - quick check: `python run_omm.py -h` for input arguments
'''


import argparse

# OpenMM imports
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd imports
from parmed import load_file
from parmed import gromacs
from parmed import unit as u

# MDTraj imports
from mdtraj import load, reporters

parser = argparse.ArgumentParser(description='OpenMM production run')
parser.add_argument('-t', '--runTime', type=float,
                    help='simulation time in nanoseconds', default=1)
parser.add_argument('-f', '--fileroot', type=str,
                    help='name of the fileroot prepared in gromacs')
parser.add_argument('-idx1', type=int,
                    help='previous index, e.g. traj0.nc', default=0)
parser.add_argument('-idx2', type=int,
                    help='current index, e.g. traj1.nc', default=1)

args = parser.parse_args()
runTime = args.runTime
fileroot = args.fileroot
idx1 = args.idx1
idx2 = args.idx2


# Load the Gromacs files
print('Loading Gromacs files...')
gromacs.GROMACS_TOPDIR = '.'
top = load_file('%s_ions.top' % fileroot,
                xyz='%s_npt.gro' % fileroot)

# Create the OpenMM system
print('Creating OpenMM System')
system = top.createSystem(nonbondedMethod=app.PME,
                          nonbondedCutoff=9.0 * u.angstroms,
                          constraints=app.HBonds,
                          hydrogenMass=1.5 * u.amu,
                          )

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinMiddleIntegrator(
    300 * u.kelvin,       # Temperature of heat bath
    1.0 / u.picoseconds,  # Friction coefficient
    4.0 * u.femtoseconds,  # Time step
)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference.
# Or do not specify the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CUDA')
# Use mixed single/double precision
properties = dict(CudaPrecision='mixed')

# Create the Simulation object
sim = app.Simulation(top.topology, system, integrator, platform, properties)

# Set the particle positions
# sim.context.setPositions(top.positions)

# Print potential energy at time=0
# state = sim.context.getState(getEnergy=True, getForces=True)
# print(state.getPotentialEnergy())

# Minimize the energy
# print('Minimizing energy')
# sim.minimizeEnergy(maxIterations=500)

print('Restart from equilibration, loading checkpoint...')
sim.loadCheckpoint('state%i.chk' % idx1)

# Set up the reporters to report energies and coordinates every 100 steps
sim.reporters.append(
    app.StateDataReporter('thermo%i.out' % idx2, 1000, step=True, potentialEnergy=True,
                          kineticEnergy=True, temperature=True, volume=True,
                          density=True, speed=True)
)

# saving the whole system usning DCDReporter from OMM
freq_save = 250000  # 250k x 4fs = 1000 ps = 1ns
sim.reporters.append(app.DCDReporter('traj%i.dcd' % idx2, freq_save))

# using NetCDFReporter from mdtraj to select a subset of the system to write into the trajectory
# freq_save2nc = 50000    # 50k x 4fs = 200 ps

pro = load('%s_npt.gro' % fileroot)
sim.reporters.append(reporters.NetCDFReporter(
    'traj%i.nc' % idx2, 0.2 * freq_save, atomSubset=pro.topology.select("not water")))


# Save the checkpoint every 1ns
sim.reporters.append(app.CheckpointReporter(
    'checkpnt%i.chk' % idx2, freq_save))

# Run production
print('Running dynamics')

sim.step(runTime * freq_save)

# Save the checkpoint at the end of the simulation
sim.saveCheckpoint('state%i.chk' % idx2)

print('Done!')
