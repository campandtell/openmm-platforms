#!/usr/bin/env python3
import openmm as mm
import openmm.app as app
from openmm.unit import *
import numpy as np

class MDSimulation:
    def __init__(self, pdb_file, params_file=None):
        """Initialize simulation parameters."""
        self.pdb = app.PDBFile(pdb_file)
        self.params = params_file
        self.system = None
        self.simulation = None
        
    def solvate_and_neutralize(self, boxPadding=1.0*nanometers, ionicStrength=0.1*molar):
        """Solvate the protein in a water box and add ions."""
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        modeller = app.Modeller(self.pdb.topology, self.pdb.positions)
        
        # Add water box
        modeller.addSolvent(
            forcefield, 
            padding=boxPadding,
            model='tip3p',
            ionicStrength=ionicStrength
        )
        
        # Update topology and positions with solvated system
        self.pdb.topology = modeller.topology
        self.pdb.positions = modeller.positions
        
        # Save the solvated system topology
        with open('solvated_system.pdb', 'w') as f:
            app.PDBFile.writeFile(self.pdb.topology, self.pdb.positions, f)
            
        print(f'System size after solvation: {modeller.topology.getNumAtoms()} atoms')
        print(f'Saved solvated system topology to solvated_system.pdb')
        
    def setup_system(self, force_field='amber14-all.xml', water_ff='amber14/tip3p.xml'):
        """Setup the system with specified force field."""
        if self.params:
            forcefield = app.ForceField(force_field, water_ff, self.params)
        else:
            forcefield = app.ForceField(force_field, water_ff)
            
        self.system = forcefield.createSystem(
            self.pdb.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*nanometers,
            constraints=app.HBonds,
            rigidWater=True
        )
        
    def setup_simulation(self, platform_name='CUDA', 
                        temperature=300*kelvin,
                        friction=1/picosecond,
                        timestep=0.002*picoseconds):
        """Setup simulation with specified platform and parameters."""
        platform = mm.Platform.getPlatformByName(platform_name)
        
        # Use Langevin Middle integrator for better temperature control
        # The friction coefficient determines how strongly the system is coupled to the heat bath
        integrator = mm.LangevinMiddleIntegrator(
            temperature,  # System temperature
            friction,     # Friction coefficient for temperature coupling
            timestep     # Integration time step
        )
        
        self.simulation = app.Simulation(
            self.pdb.topology,
            self.system,
            integrator,
            platform
        )
        self.simulation.context.setPositions(self.pdb.positions)
        
    def minimize_and_equilibrate(self):
        """Minimize and equilibrate the system."""
        print('Minimizing...')
        self.simulation.minimizeEnergy()
        
        print('Equilibrating...')
        self.simulation.context.setVelocitiesToTemperature(300*kelvin)
        self.simulation.step(10000)  # 20 ps equilibration
        
    def run_production(self, n_steps, 
                       report_interval=1000,
                       dcd_interval=1000,
                       log_file='output.txt',
                       dcd_file='trajectory.dcd'):
        """Run production simulation with specified output frequencies."""
        # Clear any existing reporters
        self.simulation.reporters.clear()
        
        # Add state data reporter for logging
        self.simulation.reporters.append(app.StateDataReporter(
            log_file, 
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=n_steps
        ))
        
        # Add trajectory reporter
        self.simulation.reporters.append(app.DCDReporter(
            dcd_file, 
            dcd_interval
        ))
        
        print(f'Running production for {n_steps} steps...')
        self.simulation.step(n_steps)
