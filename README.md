# OpenMM Simulation and Analysis Platform

A comprehensive platform for running molecular dynamics simulations using OpenMM across different computing environments: local GPU, SLURM HPC, and AWS. Includes tools for trajectory analysis, RMSF calculations, and covariance matrix generation.

## Features

- Run OpenMM simulations on:
  - Local GPU
  - SLURM HPC clusters
  - AWS GPU instances
- Automated system preparation (solvation, neutralization)
- Analysis tools:
  - CA-only trajectory extraction
  - RMSF calculation and visualization
  - Position-position covariance matrices with sliding windows

## Installation

### Local Development Setup

```bash
# Create conda environment
conda create -n openmm-dev python=3.9
conda activate openmm-dev

# Install dependencies
conda install -c conda-forge openmm mdtraj numpy matplotlib

# Clone repository
git clone https://github.com/yourusername/openmm-platforms.git
cd openmm-platforms

# Install in development mode
pip install -e .
```

### Docker Setup

```bash
# Build Docker image
docker build -t openmm-platforms .

# If you want to run GUI applications (e.g., image viewers) inside the container
# so that windows pop up on your host machine, first allow local root connections:
# (Run this on your host, not inside Docker. You'll need xhost installed on your host.)
xhost +local:root

# Run container with GPU support and X11 forwarding
docker run -it \
    --gpus all \
    -v "$(pwd):/workspace" \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --net=host \
    openmm-platforms bash
```

## Usage Example

### 1. Download Example System
```bash
# Create example system directory
mkdir -p example_systems/1btl
cd example_systems/1btl

# Download 1BTL structure
wget https://files.rcsb.org/download/1BTL.pdb
mv 1BTL.pdb 1btl.pdb
```

### 2. Run Simulation Locally

```bash
# Run simulation with local GPU
python3 scripts/local_run.py --pdb example_systems/1btl/1btl.pdb --steps 50000
```

This will:
- Load the 1BTL structure
- Solvate and neutralize the system
- Save the prepared system as 'solvated_system.pdb'
- Run 50,000 steps of simulation
- Output:
  - trajectory.dcd: Simulation trajectory
  - simulation_output.txt: Energy, temperature logs

### 3. Analyze Results

```bash
# Run analysis
python3 src/analysis.py \
    --traj trajectory.dcd \
    --top solvated_system.pdb \
    --ca-out ca_only.dcd \
    --rmsf-plot rmsf.png \
    --cov-dir covariance_results \
    --window-size 10 \
    --stride 5
```

This will generate:
- ca_only.dcd: CA-only trajectory
- rmsf.png: RMSF plot
- covariance_results/: Directory with covariance matrices

## Running on Different Platforms

### SLURM HPC

1. Edit SLURM settings in scripts/slurm_submit.sh if needed
2. Submit job:
```bash
sbatch scripts/slurm_submit.sh example_systems/1btl/1btl.pdb 500000
```

### AWS

1. Launch GPU instance:
```bash
# Launch p3.2xlarge instance with Amazon Linux 2
aws ec2 run-instances \
    --image-id ami-0123456789abcdef0 \
    --instance-type p3.2xlarge \
    --key-name your-key-pair \
    --security-group-ids sg-0123456789abcdef0
```

2. Setup instance:
```bash
# SSH into instance
ssh -i your-key.pem ec2-user@your-instance-ip

# Run setup script
bash scripts/aws_setup.sh
```

3. Run simulation:
```bash
python3 scripts/local_run.py --pdb example_systems/1btl/1btl.pdb --steps 500000
```

## Directory Structure

```
openmm-platforms/
├── README.md
├── requirements.txt
├── Dockerfile
├── src/
│   ├── __init__.py
│   ├── simulation.py
│   ├── utils.py
│   └── analysis.py
├── scripts/
│   ├── local_run.py
│   ├── slurm_submit.sh
│   └── aws_setup.sh
├── tests/
│   ├── __init__.py
│   └── test_simulation.py
└── example_systems/
    └── 1btl/
        └── 1btl.pdb
```


## License

Distributed under the MIT License. See `LICENSE` for more information.
