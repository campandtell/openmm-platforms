#!/bin/bash

# Install CUDA drivers and toolkit
sudo apt-get update
sudo apt-get install -y nvidia-driver-470 nvidia-cuda-toolkit

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/miniconda3/bin/activate

# Create and setup conda environment
conda create -n openmm-dev python=3.9 -y
conda activate openmm-dev
conda install -c conda-forge openmm mdtraj -y

# Clone repository
git clone https://github.com/yourusername/openmm-platforms.git
cd openmm-platforms
pip install -e .
