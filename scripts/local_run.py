import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.simulation import MDSimulation
from openmm.unit import *
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--params', help='Optional parameter file')
    parser.add_argument('--steps', type=int, default=500000, 
                       help='Number of simulation steps')
    args = parser.parse_args()
    
    # Initialize simulation
    sim = MDSimulation(args.pdb, args.params)
    
    # Add water and ions
    print("Adding water and ions...")
    sim.solvate_and_neutralize(
        boxPadding=1.2*nanometers,
        ionicStrength=0.15*molar
    )
    
    # Setup system and simulation
    print("Setting up system...")
    sim.setup_system()
    
    print("Setting up simulation...")
    sim.setup_simulation(
        temperature=310*kelvin,
        friction=5/picosecond,
        timestep=0.002*picoseconds
    )
    
    # Minimize and equilibrate
    print("Minimizing and equilibrating...")
    sim.minimize_and_equilibrate()
    
    # Run production
    print("Running production...")
    sim.run_production(
        n_steps=args.steps,
        report_interval=1000,
        dcd_interval=1000,
        log_file='simulation_output.txt',
        dcd_file='trajectory.dcd'
    )

if __name__ == '__main__':
    main()
