import mdtraj as md
import numpy as np
from typing import Union, Tuple, List
import matplotlib.pyplot as plt
from pathlib import Path

class TrajectoryAnalysis:
    def __init__(self, trajectory_file: str, topology_file: str):
        """
        Initialize trajectory analysis.
        
        Args:
            trajectory_file: Path to trajectory file (dcd, xtc, etc.)
            topology_file: Path to topology file (pdb, etc.)
        """
        self.traj = md.load(trajectory_file, top=topology_file)
        self.topology = md.load(topology_file)
        self.ca_traj = None
        print(f"Loaded trajectory with {self.traj.n_frames} frames and {self.traj.n_atoms} atoms")
        
    def create_ca_trajectory(self, output_file: str = None) -> md.Trajectory:
        """
        Extract CA-only trajectory.
        
        Args:
            output_file: Optional path to save reduced trajectory
            
        Returns:
            CA-only trajectory
        """
        # Select CA atoms
        ca_idx = self.traj.topology.select('name CA')
        self.ca_traj = self.traj.atom_slice(ca_idx)
        
        # Save trajectory if requested
        if output_file:
            self.ca_traj.save(output_file)
            print(f"Saved CA-only trajectory with {self.ca_traj.n_atoms} atoms to {output_file}")
        
        return self.ca_traj
    
    def calculate_rmsf(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate RMSF for CA atoms.
        
        Returns:
            Tuple of (residue indices, RMSF values)
        """
        if self.ca_traj is None:
            self.create_ca_trajectory()
            
        # Align trajectory to first frame
        aligned_traj = self.ca_traj.superpose(self.ca_traj, frame=0)
        
        # Calculate RMSF
        xyz = aligned_traj.xyz
        mean_xyz = xyz.mean(axis=0)
        diff = xyz - mean_xyz
        rmsf = np.sqrt(np.mean(np.sum(diff * diff, axis=2), axis=0))
        
        # Get residue indices
        residue_indices = [atom.residue.index for atom in self.ca_traj.topology.atoms]
        
        return np.array(residue_indices), rmsf * 10  # Convert to Angstroms
    
    def plot_rmsf(self, output_file: str = None):
        """
        Plot RMSF values for CA atoms.
        
        Args:
            output_file: Optional path to save plot
        """
        residue_indices, rmsf = self.calculate_rmsf()
        
        plt.figure(figsize=(10, 6))
        plt.plot(residue_indices, rmsf, '-o')
        plt.xlabel('Residue Index')
        plt.ylabel('RMSF (Å)')
        plt.title('Root Mean Square Fluctuation (CA atoms)')
        plt.grid(True)
        
        if output_file:
            plt.savefig(output_file)
            print(f"Saved RMSF plot to {output_file}")
        plt.show()
        
    def calculate_covariance_windows(self, 
                                   window_size: int,
                                   stride: int,
                                   output_dir: str) -> List[np.ndarray]:
        """
        Calculate CA position-position covariance matrices over sliding windows.
        
        Args:
            window_size: Number of frames per window
            stride: Number of frames to move window
            output_dir: Directory to save covariance matrices
            
        Returns:
            List of covariance matrices
        """
        if self.ca_traj is None:
            self.create_ca_trajectory()
            
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        n_frames = self.ca_traj.n_frames
        n_windows = (n_frames - window_size) // stride + 1
        
        # Calculate mass matrix once (for CA atoms)
        masses = np.array([atom.element.mass for atom in self.ca_traj.topology.atoms])
        mass_matrix = np.sqrt(np.outer(masses, masses))
        
        covariance_matrices = []
        
        for i in range(n_windows):
            start_frame = i * stride
            end_frame = start_frame + window_size
            
            print(f"Processing window {i+1}/{n_windows} (frames {start_frame}-{end_frame})")
            
            # Get coordinates for window
            window_traj = self.ca_traj[start_frame:end_frame]
            
            # Align to first frame of window
            window_traj.superpose(window_traj, 0)
            window_xyz = window_traj.xyz
            
            # Reshape to 3N coordinates
            n_atoms = window_xyz.shape[1]
            coords_3n = window_xyz.reshape(window_size, 3 * n_atoms)
            
            # Calculate covariance
            coords_mean = coords_3n.mean(axis=0)
            coords_centered = coords_3n - coords_mean
            cov_matrix = np.dot(coords_centered.T, coords_centered) / (window_size - 1)
            
            # Apply mass weighting
            mass_matrix_3n = np.kron(mass_matrix, np.ones((3, 3)))
            cov_matrix = cov_matrix * mass_matrix_3n
            
            # Save matrix
            output_file = output_dir / f'covariance_window_{start_frame}_{end_frame}.dat'
            np.savetxt(output_file, cov_matrix)
            print(f"Saved covariance matrix to {output_file}")
            
            covariance_matrices.append(cov_matrix)
            
        return covariance_matrices

def main():
    """Example usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Trajectory Analysis Tools')
    parser.add_argument('--traj', required=True, help='Trajectory file')
    parser.add_argument('--top', required=True, help='Topology file')
    parser.add_argument('--ca-out', help='Output CA-only trajectory')
    parser.add_argument('--rmsf-plot', help='Output RMSF plot file')
    parser.add_argument('--cov-dir', help='Directory for covariance matrices')
    parser.add_argument('--window-size', type=int, default=100, help='Window size for covariance')
    parser.add_argument('--stride', type=int, default=50, help='Stride for moving window')
    
    args = parser.parse_args()
    
    # Initialize analysis
    analysis = TrajectoryAnalysis(args.traj, args.top)
    
    # Create CA trajectory first
    if args.ca_out:
        analysis.create_ca_trajectory(args.ca_out)
    else:
        analysis.create_ca_trajectory()
    
    # Generate RMSF plot if requested
    if args.rmsf_plot:
        analysis.plot_rmsf(output_file=args.rmsf_plot)
    
    # Calculate covariance matrices if requested
    if args.cov_dir:
        analysis.calculate_covariance_windows(args.window_size, args.stride, args.cov_dir)

if __name__ == '__main__':
    main()
