# Use NVIDIA CUDA base image
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy \
    matplotlib \
    mdtraj \
    openmm

# Copy project files
COPY . /workspace/

# Install project in development mode
RUN pip3 install -e .

# Set default command
CMD ["/bin/bash"]
